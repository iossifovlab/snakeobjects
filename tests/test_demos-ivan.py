import os,glob
import pytest

DEMOS_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    'demos',
    )

ALL_DEMOS = list(glob.glob(DEMOS_DIR + "/*"))

# print("DEMOS_DIR", DEMOS_DIR)
# print("ALL_DEMOS", ALL_DEMOS)

@pytest.mark.parametrize("DD", ALL_DEMOS)
def test_bla(DD,tmpdir):
    print("DDDDDDD",DD,tmpdir)
    assert os.system(f'cp -r {DD}/* {tmpdir}') == 0
    assert os.system(f'ls -l {tmpdir}') == 0
    assert os.system(f'{tmpdir}/run.sh') == 0

