from pybustools.pybustools import iterate_bus_cells_umi_multiple
from pybustools.pybustools import Bus
from pybustools.parallel_generators import ParallelCellUMIGenerator
import collections
import toolz
import numpy as np
import pandas as pd
import tqdm
import pickle
from phantom import emit_records_based_on_gene, phantom_create_dataframes, phantom_prep_binom_regr


BUSDIR = '/home/mstrasse/mountSSD/kallisto_out'
# flowcell V300026370
# flowcell = 'dnbseqg400.V300026370'
# lanes = {'L1': [
#     'dnbseqg400.V300026370_88A.E11B_2-718745_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
#     'dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
#     'dnbseqg400.V300026370_88A.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
# ]}

# flowcell V300039753
flowcell = 'dnbseqg400.V300039753'
lanes = {
    'L1': ['dnbseqg400.V300039753.E15A_2-750715_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs', 'dnbseqg400.V300039753.E15B_2-750716_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs'],
    'L2': ['dnbseqg400.V300039753.E15C_2-751611_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs', 'dnbseqg400.V300039753.E15D_2-751118_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs'],
    'L3': ['dnbseqg400.V300039753.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs', 'dnbseqg400.V300039753.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs'],
    }


# # novaseq flowcell
# flowcell = 'novaseq.190919_A00266_0278_BHFJY5DRXX'
# lanes = {'L1': ['novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs',
#                 'novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs']
#         }

BUSDIR = '/home/mstrasse/mount/kallisto_out'
flowcell = 'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC'
lanes = {
'L1': [
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E07B_de-cryo_2-718742_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E08A_de-cryo_2-718743_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E08C_de-cryo_2-718744_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E10A_2-718739_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E10B_2-718740_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E10D_2-718741_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E11B_2-718745_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E11C_2-718746_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E11D_2-718747_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E11E_2-718748_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E12B_2-729661_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E12C_2-729662_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.E12D_2-729663_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.L06C_2-718733_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.L07A_2-718734_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.L07B_2-718735_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.L07C_2-718736_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.L07D1_2-718737_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.191218_A00266_0321_BHYJCLDSXX_10XmixteGrundbergATAC.L07D2_2-718738_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
]
}

lane = 'L1'
samples = lanes[lane]
# busfiles = [f'{BUSDIR}/{s}/kallisto/sort_bus/bus_output/output.corrected.sort.bus' for s in samples]
bus = {s: Bus(f'{BUSDIR}/{s}/kallisto/sort_bus/bus_output/') for s in samples}

# quick check: are all these EC-dicts compatible? they are not equal for sure!
# EC = {}
# for s, b in bus.items():
#     for ec, transcript_list in b.ec_dict.items():
#         if not ec in EC:
#             EC[ec] = transcript_list
#         else:
#             assert EC[ec] == transcript_list

df_finger, conditional = phantom_create_dataframes(bus)
with open(f'/tmp/{flowcell}-{lane}_phantomPurger.pkl', 'wb') as fh:
    pickle.dump([df_finger, conditional], fh)


"""
=============================================================================
create the dataframes for the purger
This is the comp. intensive part
=============================================================================
"""
for lane, samples in lanes.items():
    bus = {s: Bus(f'{BUSDIR}/kallisto_out/{s}/kallisto/sort_bus/bus_output/') for s in samples}

    df_finger, conditional = phantom_create_dataframes(bus)
    with open(f'/tmp/{flowcell}-{lane}_phantomPurger.pkl', 'wb') as fh:
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

with open(f'/tmp/{flowcell}_chimer_cells.pkl', 'wb') as fh:
    pickle.dump([chimer_per_cell, chimer_per_cell_notorigin], fh)

# ==========================================================
# ==========================================================
# ==========================================================
# count the UMI overlap of cells
# ===========================================================
samples = [
    'dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300026370_88A.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300039753.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'dnbseqg400.V300039753.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs',
    'novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs',
    'novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs'
]
busfiles = [f'{BUSDIR}/kallisto_out/{s}/kallisto/sort_bus/bus_output/output.corrected.sort.bus' for s in samples]
bus_dict = {n:b for n,b in zip(samples, busfiles)}
import itertools
import tqdm
from pybustools.parallel_generators import ParallelCellGenerator
from pybustools.pybustools import iterate_bus_cells_multiple

if True:  # run parallel
    PG = ParallelCellGenerator(bus_dict, decode_seq=False, queue_size=10000)
    PG.start_queues()
    I = PG.iterate()
else:  # run serial
    I = iterate_bus_cells_multiple(samples, busfiles, decode_seq=False)


def jacard(list1, list2):
    s1 = set(list1)
    s2 = set(list2)
    return len(s1 & s2) / len(s1 | s2)

cb_jaccard = []
counter = 0
for cb, info_dict in tqdm.tqdm(I):
    jac = {}
    for p1, p2 in itertools.combinations(samples, 2):
        if p1 in info_dict and p2 in info_dict:
            umi1 = [_[0] for _ in info_dict[p1]]
            umi2 = [_[0] for _ in info_dict[p2]]
            # umi2 = info_dict[p2]
            j = jacard(umi1, umi2)
            jac[(p1,p2)] = j
    jac['CB'] = cb
    cb_jaccard.append(jac)
    counter+=1
# parallel: 2504888it [19:41, 2120.96it/s]
# serial :  2504888it [24:08, 1729.65it/s]
