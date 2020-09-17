import pytest, os, sys
from glob import glob

cwd = os.getcwd()


def demo_test(d):
    os.chdir('demo/'+d)
    os.system('./clean.sh')
    os.system('./run.sh')

    os.chdir('test-out')
    fnames = glob('*/*/*')
    os.chdir('..')
    D = ' objLinks/'
    for fn in fnames:
        assert os.system('diff test-out/' + fn + D + fn) == 0
    os.chdir('../..')

def test_d1():
    demo_test('d1')

def test_d3():
    demo_test('d3')

