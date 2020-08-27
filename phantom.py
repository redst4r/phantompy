from pybustools.pybustools import iterate_bus_cells_umi_multiple
from pybustools.pybustools import Bus
from pybustools.parallel_generators import ParallelCellUMIGenerator
import collections
import toolz
import numpy as np
import pandas as pd
import tqdm


# def info_list_to_genes(info_list, samplename):
#     """
#     each iteration return (cb, umi), list(tuples)
#
#     usually the list is len(1): i.e. a molecule mapping to a EC.
#     However, sometimes theres multiple reads from the same molecule that map to different ECs
#     """
#     genes = [set(bus[samplename].ec_dict[ec]) for ec, count, flag in info_list]
#     # since its a single molecule, there MUST be some overlap in the genes
#     o = genes[0]
#     for s in genes[1:]:
#         o = o & s
#     return o
TMPbus = collections.namedtuple("TMPbus", 'cb umi ec counts flag samplename')


class DisjointSets():
    """
    d = DisjointSets()
    d.add_set({1,2,3}, 'A')
    d.add_set({4,5}, 'B')
    d.add_set({4,7}, 'X')  # this will merge with B
    d.add_set({1,7}, 'Y')  # this should merge the entire thing
    """

    def __init__(self):
        self.disjoint_sets = {}

    def add_set(self, aset, name):

        # look for disjoint sets that share a member with the current set
        candidate_setix = []
        for n, s in self.disjoint_sets.items():
            if len(s & aset) > 0:  # there;s a shared element
                candidate_setix.append(n)

        # now, anything thats in the candidate set (+aset iself)
        # for a new aggregated disjoint set!
        newset = aset
        # new_name = name
        new_name = (name,)  # turn into tuple
        for n in candidate_setix:
            newset = newset | self.disjoint_sets[n]
            # new_name = new_name + '_' + n
            new_name = new_name + n
            del self.disjoint_sets[n]

        self.disjoint_sets[new_name] = newset

    def n_disjoint_sets(self):
        return len(self.disjoint_sets)


def emit_records_based_on_gene(cb, umi, info_dict, busobject_dict):
    """
    this splits the record into groups that map to the same gene/transcript

    info_dict contains muliple bus entries from different samples.
    we want to group theese such that bus-records coming form the same
    molecules are groupd
    """
    bus_list = []

    for samplename, info_list in info_dict.items():
        for entry in info_list:
            _t = TMPbus(cb, umi, entry[0], entry[1], entry[2], samplename)
            bus_list.append(_t)

    # shortcut: if its just a single bus entry, yield the original cb,umi,dict
    if len(bus_list) == 1:
        yield (cb, umi), info_dict
    else:
        # get the genes corresponding to each record
        genes = {}
        for b in bus_list:
            bus_object = busobject_dict[b.samplename]
            # TODO this ASSUMES THAT all transcript IDs encode the same stuff across busfiles!!
            g = bus_object.ec_dict[b.ec]
            genes[b] = set(g)

        # group the samples, based on if they have overlaping genes
        DS = DisjointSets()
        for b, gene_set in genes.items():
            DS.add_set(gene_set, name=b)

        # now, each dijoint set is identified by a tuple of samples
        for bus_tuple in DS.disjoint_sets.keys():
            # the keys are a consistent single molecule
            # so we must create only ONE result per key
            # if they are from the same sample, add up reads
            emitted_dict = {}
            for b in bus_tuple:
                if b.samplename in emitted_dict:
                    ec, counts, flag = emitted_dict[b.samplename]
                    emitted_dict[b.samplename] = (ec, counts+b.counts, flag)
                else:
                    fake_ec = -1
                    emitted_dict[b.samplename] = (fake_ec, b.counts, b.flag)
            emitted_dict = toolz.valmap(lambda x: [x], emitted_dict)  # expected to yield a list
            yield (cb, umi), emitted_dict


def _bus_check_transcript_equivalence(list_of_busobjects):
    # quick check: are all these transcript-dicts compatible?
    # WE ASSUME THIS IN `emit_records_based_on_gene`
    for bus in list_of_busobjects:
        assert bus.transcript_dict == list_of_busobjects[0].transcript_dict, \
            "transcript ids map to different transcripts!!"


def _create_fingerprint(info_dict):
    """
    turns an info_dict (a dict of bus-records from different experiments all
    corresponding to the same molecule) into a fingerprint:
    counting the #reads per experiment of that molecule
    """
    fingerprint = collections.defaultdict(int)
    for name, info_list in info_dict.items():
        # for each sample, count the number of reads
        n_reads = 0
        for ec, counts, flag in info_list:
            n_reads += counts
        fingerprint[name] = n_reads
    return fingerprint


def _counter_to_df(fingerprints_counter, samples):
    k, v = zip(*fingerprints_counter.items())
    df_finger = pd.DataFrame(k, columns=samples)
    k_chimera = np.sum(df_finger.values > 0, 1)
    r = np.sum(df_finger.values, 1)

    df_finger['k_chimera'] = k_chimera
    df_finger['r'] = r
    df_finger['freq'] = v

    # calc the conditional (on r) prob of finding a molecule
    # ampliefied r times in each experiment
    conditional = []

    for r in df_finger['r'].unique():
        _tmp = df_finger.query('r==@r')

        # each sample has k amount of reads.
        # but we also obersve the entire fingerprint n times
        # so each expeiment really had k*n amount of reads
        v = {}
        for s in samples:
            v[s] = np.sum(_tmp[s] * _tmp['freq'])
        v_total = sum(v.values())
        v = toolz.valmap(lambda x: x/v_total, v)
        for s in samples:
            conditional.append({'r': r, 'sample': s, 'fraction': v[s]})
    conditional = pd.DataFrame(conditional).pivot_table(index='r', columns='sample', values='fraction')

    return df_finger, conditional


def phantom_create_dataframes(busobject_dict):
    """
    main function of PhantomPurger: for the set of samples (dict of busfiles),
    determine the molecules appearing in multiple experiments.
    Instead of trackign every molecule, we just do bookkeeping on the
    "experiment" fingerprints of each molecule:
     - each molecule gets a charateristic fingerprint/vector
       (hhow often it occurs in which experiment)
     - instead of creating one entry per molecule, we just keep track of the
       frequencies of these fingerprints
     - With four samples, a molecule might have the fingerprint: [0, 10, 2, 1]
       (doesnt occur in exp1, 10x in exp2, 2x in exp3...)
    """

    # unpack the dict; its easier
    samples, busfiles = [], []
    for s, bus in busobject_dict.items():
        samples.append(s)
        busfiles.append(str(bus.bus_file))

    _bus_check_transcript_equivalence(list(busobject_dict.values()))

    # initialie iterators
    PARALLEL = False
    if PARALLEL:
        PG = ParallelCellUMIGenerator({s: b for s, b in zip(samples, busfiles)},
                                      decode_seq=False,
                                      queue_size=10000)
        PG.start_queues()
        bus_iter = PG.iterate()
    else:
        bus_iter = iterate_bus_cells_umi_multiple(samples, busfiles, decode_seq=False)

    # keep track of how molecules get amplified in the different samples
    amp_factors_per_sample = {n: collections.defaultdict(int) for n in samples}
    fingerprints_counter = collections.defaultdict(int)

    for (cb_, umi_), info_dict_ in tqdm.tqdm(bus_iter):
        # THIS SPLIts according to same gene
        for (cb, umi), info_dict in emit_records_based_on_gene(cb_, umi_, info_dict_, busobject_dict):
            fingerprint = _create_fingerprint(info_dict)
            for name, n_reads in fingerprint.items():
                amp_factors_per_sample[name][n_reads] += 1

            fp = [fingerprint[_] for _ in samples]
            fingerprints_counter[tuple(fp)] += 1

    # if PARALLEL:
    #     PG.cleanup()

    df_finger, conditional = _counter_to_df(fingerprints_counter, samples)
    return df_finger, conditional


def _group_by_ec(cb, umi, info_dict):

    bus_list = []
    for samplename, info_list in info_dict.items():
        for entry in info_list:
            _t = TMPbus(cb, umi, entry[0], entry[1], entry[2], samplename)
            bus_list.append(_t)

    groups = collections.defaultdict(list)
    for r in bus_list:
        groups[r.ec].append(r)

    for ec, record_list in groups.items():
        emit_dict = {}
        for r in record_list:
            if r.samplename in emit_dict:
                emit_dict[r.samplename].append((r.ec, r.counts, r.flag))
            else:
                emit_dict[r.samplename] = [(r.ec, r.counts, r.flag)]
        yield (cb, umi), emit_dict


def phantom_prep_binom_regr(df_finger):
    """
    turn the fingerprint dataframe into a DF that we can use in binomial
    regression. Summarises the number of chimeric/nonchimeric molecule as a
    function of r (r=sum of reads of the molecule)
    """
    # number of non-chimeric molecules as a function of r
    df_non_chimers = df_finger.query('k_chimera==1')
    df_non_chimers = df_non_chimers.groupby('r')[['freq']].sum()
    df_non_chimers.rename({'freq': 'freq_non_chimer'}, axis=1, inplace=True)
    df_non_chimers.head()

    # TODO not sure if the k==2 shouldnt be k>=2!!!!
    import warnings
    warnings.warn('TODO not sure if the k==2 shouldnt be k>=2!!!!')
    df_chimers = df_finger.query('k_chimera==2').groupby('r')[['freq']].sum()
    df_chimers.rename({'freq': 'freq_chimer'}, axis=1, inplace=True)
    df_chimers.head()

    df_binom = df_non_chimers.merge(df_chimers, left_index=True, right_index=True, how='outer')
    df_binom = df_binom.replace({np.nan: 0})

    return df_binom
