import pandas as pd
import numpy as np
import sys
import fire
import pickle
import tqdm
from pybustools.pybustools import iterate_bus_cells_umi_multiple
from pybustools.busio import write_busfile, get_header_info, read_binary_bus
from pybustools.utils import h5_to_bus
import collections
VERBOSE = True


def h5_convert(h5file, busfile, TMPDIR):
    """
    convert cellranger molecule.h5 files to the bus format.
    This is used in phantompurger for those data where not all samples are available as bam/fastq files
    
    """
    h5_to_bus(h5file, busfile, TMPDIR)


def identify_suspicious(outfile, *busfiles):
    """
    Identify suspicious CUGs (that probably swapped). Done by looking at all
    datasets from a singl run/lane and identify CUGs that occur more than once
    """
    if VERBOSE:
        print(f'Processing {len(busfiles)} busfiles:')
        for b in busfiles:
            print(f'\t{b}')

    samples = busfiles  # just dummy names for each busfile
    bus_iter = iterate_bus_cells_umi_multiple(samples, busfiles,
                                              decode_seq=False)

    dubious_cb_umi = []

    # also keep track of how much the experiments overlapped in terms of umis
    # key are all potential overlaps (for 3 experiments: 12, 13, 23, 123)
    overlaps = collections.defaultdict(int)

    for (cb, umi), info_dict_ in tqdm.tqdm(bus_iter):
        if len(info_dict_) > 1:
            # this molecule occured in more then 1 samples,
            # hence a hopped candidate
            dubious_cb_umi.append((cb, umi))
        s = frozenset(info_dict_.keys())
        overlaps[s] += 1

    with open(outfile, 'wb') as fh:
        pickle.dump(dubious_cb_umi, fh)

    print(f'{len(dubious_cb_umi)} CUGs marked suspicious')

    print('Overlaps')
    for k, v in overlaps.items():
        print(f'{k}\t\t\t\t: {v}')


def filter_busfile(inbus, outbus, suspicious):
    """
    filtering the inbus file for CUGs that are listed in the suspicious-file
    """
    if VERBOSE:
        print(f'Filtering {inbus} using {suspicious}')

    with open(suspicious, 'rb') as fh:
        dubious_cb_umi = pickle.load(fh)
        dubious_cb_umi = set(dubious_cb_umi)

    _, cb_len, umi_len, _ = get_header_info(inbus)

    def _gen():
        n_filtered = 0
        n_total = 0
        for record in tqdm.tqdm(read_binary_bus(inbus, decode_seq=False)):
            n_total += 1
            if (record.CB, record.UMI) in dubious_cb_umi:
                n_filtered += 1
                continue
            yield record

        print(f'{n_filtered}/{n_total} ({100 * n_filtered/n_total:.3f}%) CUGs filtered')

    write_busfile(outbus, _gen(), cb_len, umi_len)


if __name__ == '__main__':
    fire.Fire({
        'suspicious': identify_suspicious,
        'filter': filter_busfile,
        'h5convert': h5_convert
    })


"""
python cruk-phantom-cli.py suspicious  /tmp/out.pkl \
  /home/mstrasse/TB4/kallisto_out_trimmed/dnbseqg400.V300039753.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs/kallisto/sort_bus/bus_output/output.corrected.sort.bus \
  /home/mstrasse/TB4/kallisto_out_trimmed/dnbseqg400.V300039753.L05B_2-658953_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs/kallisto/sort_bus/bus_output/output.corrected.sort.bus
"""
