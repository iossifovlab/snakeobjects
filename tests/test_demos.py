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
    os.system('./clean.sh')
    os.chdir('../..')

def test_d1():
    demo_test('d1')

def test_d2():
    os.chdir('demo/d2')
    os.system('./clean.sh')
    os.system('./run.sh')
    for prj in ['projA', 'projB']:
        os.chdir(prj)
        os.chdir('test-out')
        fnames = glob('*/*/*')
        os.chdir('..')
        D = ' objLinks/'
        for fn in fnames:
            assert os.system('diff test-out/' + fn + D + fn) == 0
        os.chdir('..')
    os.system('./clean.sh')
    os.chdir('../..')

def test_d3():
    demo_test('d3')

def test_d4():
    demo_test('d4')

def test_d5():
    demo_test('d5')

def test_d6():
    demo_test('d6')

def test_d7():
    os.chdir('demo/d7')
    os.system('./clean.sh')
    os.system('./run.sh')
    assert os.system('diff correct_test_out.txt test_out.txt') == 0
    os.system('./clean.sh')
    os.chdir('../..')

def test_d8():
    demo_test('d8')


