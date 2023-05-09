from pathlib import Path

from dpf_mrcfile import MRC


def test_load_mrc_01(test_data_dir: Path):
    test_mrc = test_data_dir / "EMD-3001.map"
    m = MRC()
    m.read_mrc(test_mrc)
    assert len(m.volume_data) != 0
