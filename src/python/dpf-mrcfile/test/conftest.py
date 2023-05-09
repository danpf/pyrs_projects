from distutils import dir_util
import os

import pytest


@pytest.fixture
def test_data_dir(tmpdir, request):
    root_path = request.module.__file__
    root_dir = os.path.dirname(root_path)
    input_dir = "{}/{}".format(root_dir, "test_data")
    dir_util.copy_tree(input_dir, tmpdir.strpath)
    return tmpdir
