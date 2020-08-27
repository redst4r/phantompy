import pytest
from pybustools import busio
from pybustools.pybustools import iterate_bus_cells_umi_multiple, Bus
from phantom import emit_records_based_on_gene, phantom_create_dataframes
import os


def test_emit_records(tmp_path, ec_matrix_file, transcript_file):
    """
    two busfiles, they have the same CB/UMI, but it maps to different EC/genes
    emit_records_based_on_gene should yield multiple entries!!
    """
    cb_length = 4
    umi_length = 3

    records1 = [  # #      CB     UMI    EC COUNT FLAG
        busio.Bus_record('ATAT', 'GGG', 1, 10, 1),
    ]
    fname1 = tmp_path / 'some1.bus'
    busio.write_busfile(fname1, records1, cb_length, umi_length)
    bus1 = Bus(folder='/', bus_name=fname1, ec_name=ec_matrix_file, transcript_name=transcript_file)

    records2 = [
        busio.Bus_record('ATAT', 'GGG', 9, 20, 1),
    ]
    fname2 = tmp_path / 'some2.bus'
    busio.write_busfile(fname2, records2, cb_length, umi_length)
    bus2 = Bus(folder='/', bus_name=fname2, ec_name=ec_matrix_file, transcript_name=transcript_file)

    busobject_dict = {'s1': bus1, 's2': bus2}
    bus_iter = iterate_bus_cells_umi_multiple(['s1', 's2'], [fname1,  fname2], decode_seq=False)
    for (cb_, umi_), info_dict_ in bus_iter:
        # THIS SPLIts according to same gene
        print('----------')
        print(info_dict_)
        print('----------')

        counter = 0
        for (cb, umi), info_dict in emit_records_based_on_gene(cb_, umi_, info_dict_, busobject_dict):
            counter += 1
            print(cb, umi)
            print(info_dict)
        assert counter == 2


def test_fingerprint(tmp_path, ec_matrix_file, transcript_file):
    """
    This should create 3 fingerprints:
    [10, 0]  x 2
    [99, 0]  x 1
    [0, 20]  x 1
    """
    cb_length = 4
    umi_length = 3

    records1 = [  # #      CB     UMI    EC COUNT FLAG
        busio.Bus_record('ATAT', 'GGG', 1, 10, 1),
        busio.Bus_record('ATAT', 'TTT', 1, 10, 1),
        busio.Bus_record('CTAT', 'GGG', 1, 99, 1),
    ]
    fname1 = tmp_path / 'some1.bus'
    busio.write_busfile(fname1, records1, cb_length, umi_length)
    bus1 = Bus(folder='/', bus_name=fname1, ec_name=ec_matrix_file, transcript_name=transcript_file)

    records2 = [
        busio.Bus_record('ATAT', 'GGG', 9, 20, 1),
    ]
    fname2 = tmp_path / 'some2.bus'
    busio.write_busfile(fname2, records2, cb_length, umi_length)
    bus2 = Bus(folder='/', bus_name=fname2, ec_name=ec_matrix_file, transcript_name=transcript_file)

    busobject_dict = {'s1': bus1, 's2': bus2}

    df, _cond = phantom_create_dataframes(busobject_dict)

    print(df)
    assert len(df) == 3
    assert df['freq'].sum() == 4


@pytest.fixture
def ec_matrix_file():
    "creates an ec_file with 10 entries"
    import tempfile
    fname = tempfile.mktemp()
    with open(fname, 'w') as fh:
        for i in range(10):
            fh.write(f'{i} {i}\n')
    yield fname
    os.remove(fname)


@pytest.fixture
def transcript_file():
    "creates an transcript_file with 10 entries"
    import tempfile
    fname = tempfile.mktemp()
    with open(fname, 'w') as fh:
        for i in range(10):
            fh.write(f'ENST00000000000{i}.1\n')
    yield fname
    os.remove(fname)
