from pybustools.pybustools import iterate_bus_cells_umi_multiple
import collections
import toolz
from pybustools.pybustools import Bus
import numpy as np
import pandas as pd
import tqdm
import pickle


def info_list_to_genes(info_list, samplename):
    """
    each iteration return (cb, umi), list(tuples)

    usually the list is len(1): i.e. a molecule mapping to a EC.
    However, sometimes theres multiple reads from the same molecule that map to different ECs
    """
    genes = [set(bus[samplename].ec_dict[ec]) for ec, count, flag in info_list]
    # since its a single molecule, there MUST be some overlap in the genes
    o = genes[0]
    for s in genes[1:]:
        o = o & s
    return o


class DisjointSets():

    def __init__(self):
        self.disjoint_sets = {}

    def add_set(self, aset, name):

        # look for disjoint sets that share a member with the current set
        candidate_setix = []
        for n, s in self.disjoint_sets.items():
            if len(s & aset) > 0: # there;s a shared element
                candidate_setix.append(n)

        # now, anything thats in the candidate set (+aset iself)
        # for a new aggregated disjoint set!
        newset = aset
        # new_name = name
        new_name = (name,)  # turn into tuple
        for n in candidate_setix:
            newset = newset | self.disjoint_sets[n]
            # new_name = new_name + '_' + n
            new_name = new_name +  n
            del self.disjoint_sets[n]

        self.disjoint_sets[new_name] = newset

    def n_disjoint_sets(self):
        return len(self.disjoint_sets)
"""
d = DisjointSets()
d.add_set({1,2,3}, 'A')
d.add_set({4,5}, 'B')
d.add_set({4,7}, 'X')  # this will merge with B
d.add_set({1,7}, 'Y')  # this should merge the entire thing
"""

TMPbus = collections.namedtuple("TMPbus", 'cb umi ec counts flag samplename')


def emit_records_based_on_gene(cb, umi, info_dict):
    """
    this splits the record into groups that map to the same gene

    info_dict contains muliple bus entries from different samples.
    we want to group theese such that bus-records coming form the same molecules
    are groupd
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
            g = bus[b.samplename].ec_dict[b.ec]
            genes[b] = set(g)

        # group the samples, based on if they have overlaping genes
        DS = DisjointSets()
        for b, gene_set in genes.items():
            DS.add_set(gene_set, name=b)

        disjoint_sets = DS.disjoint_sets

        for bus_tuple in disjoint_sets.keys():
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
            emitted_dict = toolz.valmap(lambda x: [x], emitted_dict) # expected to yield a list
            yield (cb, umi), emitted_dict


def phantom_create_dataframes(samples, busfiles):

    I = iterate_bus_cells_umi_multiple(samples, busfiles, decode_seq=False)

    # keep track of how molecules get amplified in the different samples
    amp_factors_per_sample = {n: collections.defaultdict(int) for n in samples}

    """
    to keep track of the molecules:
     - each molecule gets a charateristic fingerprint/vector (hhow often it occurs in which experiment)
     - instead of creating one entry per molecule, we just keep track of the frequencies of these fingerprints
    """
    fingerprints_counter = collections.defaultdict(int)

    for (cb_, umi_), info_dict_ in tqdm.tqdm(I):
        # THIS SPLIts according to same gene
        for (cb, umi), info_dict in emit_records_based_on_gene(cb_, umi_, info_dict_):
            fingerprint = collections.defaultdict(int)
            for name, info_list in info_dict.items():
                # for each sample, count the number of reads
                n_reads = 0
                for ec, counts, flag in info_list:
                    n_reads += counts
                fingerprint[name] = n_reads
                amp_factors_per_sample[name][n_reads] += 1
        fp = [fingerprint[_] for _ in samples]
        fingerprints_counter[tuple(fp)] += 1

    k, v = zip(*fingerprints_counter.items())
    df_finger = pd.DataFrame(k, columns=samples)
    k_chimera = np.sum(df_finger.values > 0, 1)
    r = np.sum(df_finger.values, 1)

    df_finger['k_chimera'] = k_chimera
    df_finger['r'] = r
    df_finger['freq'] = v

    # calc the conditional (on r) prob of finding a molecule ampliefied r times in each experiment
    conditional = []

    for r in df_finger['r'].unique():
        _tmp = df_finger.query('r==@r')

        # each sample has k amount of reads. but we also obersve the entire fingerprint n times
        # so each expeiment really had k*n amount of reads
        v = {}
        for s in samples:
            v[s] = np.sum(_tmp[s] * _tmp['freq'])
        v_total = sum(v.values())
        v = toolz.valmap(lambda x: x/v_total, v)
        _d = {'r': r}
        for s in samples:
            conditional.append({'r': r, 'sample': s, 'fraction': v[s]})
    conditional = pd.DataFrame(conditional).pivot_table(index='r', columns='sample', values='fraction')

    return df_finger, conditional


# flowcell V300026370
flowcell = 'dnbseqg400.V300026370'
samples = [
    'dnbseqg400.V300026370_88A.E11B_2-718745_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300026370_88A.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
]

# flowcell V300039753
flowcell = 'dnbseqg400.V300039753'
samples = [
    'dnbseqg400.V300039753.E15A_2-750715_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300039753.E15B_2-750716_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300039753.E15C_2-751611_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300039753.E15D_2-751118_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300039753.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300039753.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
]

# novaseq flowcell
flowcell = 'novaseq.190919_A00266_0278_BHFJY5DRXX'
samples = [
    'novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs',
    'novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs'
]

busfiles = [f'/home/mstrasse/mountSSD/kallisto_out/{s}/kallisto/sort_bus/bus_output/output.corrected.sort.bus' for s in samples]
bus = {s: Bus(f'/home/mstrasse/mountSSD/kallisto_out/{s}/kallisto/sort_bus/bus_output/') for s in samples}

"""
=============================================================================
create the dataframes for the purger
This is the comp. intensive part
=============================================================================
"""

df_finger, conditional = phantom_create_dataframes(samples, busfiles)
with open(f'/tmp/{flowcell}_phantomPurger.pkl', 'wb') as fh:
    pickle.dump([df_finger, conditional], fh)


"""
-------------------------------------------------------------------
for each CB, flag the number of potential dubious molecules
-------------------------------------------------------------------
"""
chimer_per_cell = {s: collections.defaultdict(int) for s in samples} # this coutns the chimeric molecules, but not the potential ORIGIN of the molecule
chimer_per_cell_notorigin = {s: collections.defaultdict(int) for s in samples}  # just count the molecules that prob NOT originated in this cell

I = iterate_bus_cells_umi_multiple(samples, busfiles, decode_seq=False)
# I = toolz.take(10_000_000, I)
for (cb_, umi_), info_dict_ in tqdm.tqdm(I):

    # THIS SPLIts according to same gene
    for (cb, umi), info_dict in emit_records_based_on_gene(cb_, umi_, info_dict_):
        fingerprint = collections.defaultdict(int)
        argmax_sample = None
        argmax_read = -1
        for name, info_list in info_dict.items():
            # for each sample, count the number of reads
            n_reads = 0
            for ec, counts, flag in info_list:
                n_reads += counts
            fingerprint[name] = n_reads
            if n_reads > argmax_read:
                argmax_sample = name
                argmax_read = n_reads

        k_chim = sum([fingerprint[_] > 0 for _ in info_dict.keys()])

        if k_chim > 1:
            for s in info_dict.keys():
                chimer_per_cell[s][cb] += 1
                if s != argmax_sample:
                    chimer_per_cell_notorigin[s][cb] +=1


with open('/tmp/dnbseqg400.V300039753_chimer_cells.pkl', 'wb') as fh:
    pickle.dump([chimer_per_cell, chimer_per_cell_notorigin], fh)

# samples_to_index = {s: i for i, s in enumerate(samples)}
#
# row_ix = []
# col_ix = []
# data = []
# counter = 0
# # a table to keep track of the non-chimers (unique to experiment) as a function of the amp.factor r
# non_chimers = collections.defaultdict(int)
# chimers = collections.defaultdict(int)

# for (cb_, umi_), info_dict in tqdm.tqdm(I):
#
#     # THIS SPLIts according to same gene
#     for cb, umi in emit_records_based_on_gene(cb_, umi_, info_dict):
#         # only found in on experiment
#         if len(info_dict) == 1:
#             #however, this molecule migt still be present in reads that map to different EC
#             (name, info_list), = info_dict.items()
#
#             # genes = info_list_to_genes(info_list, name)
#             # if len(genes) == 0:
#             #     print('weird, single molecule mapping to two different genes')
#             #     continue
#
#             # this single molecule has multipel reads to it, which might also
#             # be distirbuted over several bus entries
#             n_reads = 0
#             for ec, counts, flag in info_list:
#                 n_reads += counts
#             non_chimers[n_reads] += 1
#
#             # row_ix.append(counter)
#             # col_ix.append(samples_to_index[name])
#             # data.append(n_reads)
#             amp_factors_per_sample[name][n_reads] += 1
#
#         if len(info_dict) > 1:
#             # this is a potential phantom molecule.
#             # same cb/umi and maps to the same gene across multiple experiments
#             reads_per_sample = {}
#             for name, info_list in info_dict.items():
#                 # for each sample, count the number of reads
#                 n_reads = 0
#                 for ec, counts, flag in info_list:
#                     n_reads += counts
#                 reads_per_sample[name] = n_reads
#
#                 amp_factors_per_sample[name][n_reads] += 1
#
#             n_reads_total = np.sum(list(reads_per_sample.values()))
#             chimers[n_reads_total] += 1
#
#
#             # the_row = np.zeros(len(samples))
#             # for name, reads in reads_per_sample.items():
#             #     # the_row[samples_to_index[name]] = reads
#             #     row_ix.append(counter)
#             #     col_ix.append(samples_to_index[name])
#             #     data.append(reads)
#         counter += 1
