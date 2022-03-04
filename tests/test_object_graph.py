import pytest
import os
from snakeobjects.ObjectGraph import ObjectGraph, ObjectGraphException, load_object_graph


@pytest.fixture
def OG():
    return ObjectGraph()


def test_create(OG):
    assert OG.get_object_types() == []


def test_add(OG):
    o = OG.add('t', 'o')
    assert OG.get_object_types() == ['t']
    assert OG['t', 'o'] == o
    assert OG['t', 'o'].oType == 't'
    assert OG['t', 'o'].oId == 'o'
    assert OG['t', 'o'].deps == []
    assert OG['t', 'o'].params == {}
    assert OG['t'] == [o]


def test_dep_par(OG):
    OG.add('A', '2')
    OG.add('A', '1')
    OG.add('A', '3')
    OG.add('B', 'o', deps=OG['A'])
    assert OG['B', 'o'].deps == [OG['A', '2'], OG['A', '1'], OG['A', '3']]


def test_dep_par(OG):
    OG.add('t', 'o', params={'p': 'v'})
    assert OG['t', 'o'].params == {'p': 'v'}


def test_dep_par(OG):
    OG.add('t', 'o')
    with pytest.raises(ObjectGraphException) as excinfo:
        OG.add('t', 'o')
    assert "is already added" in str(excinfo.value)


def test_otype_order(OG):
    from random import random
    inOrder = [str(random()) for _ in range(100)]
    for t in inOrder:
        OG.add(t, 'o')
    assert OG.get_object_types() == inOrder


def test_save_load(tmp_path):
    OG = ObjectGraph()
    a = OG.add("A", "o", {"gosho": "pesho", "opaa": "3"})
    b = OG.add("B", "1", deps=OG['A'])
    b = OG.add("B", "2", deps=OG['A'])
    b = OG.add("B", "3", deps=OG['A'])
    c = OG.add("C", "o", deps=OG['B'])
    OG.save(tmp_path / 'a.OG')
    OG1 = load_object_graph(tmp_path / 'a.OG')
    OG1.save(tmp_path / 'b.OG')
    assert os.system('diff %s %s' % (str(tmp_path / 'a.OG'), str(tmp_path / 'b.OG'))) == 0


def test_deepDeps_bug(OG):
    def p(name, ol):
        return print(name, ",".join([f"{o.oType}:{o.oId}" for o in ol]))

    def s(ol):
        return ",".join([f"{o.oType}:{o.oId}" for o in ol])

    OG = ObjectGraph()
    a = OG.add("A", "o")
    b = OG.add("B", "o", deps=OG['A'])
    c = OG.add("C", "o", deps=[OG['A', 'o'], OG['B', 'o']])

    assert s(c.deepDeps('A', level=2)) == "A:o"


def test_deepDeps(OG):
    def p(name, ol):
        return print(name, ",".join([f"{o.oType}:{o.oId}" for o in ol]))

    def s(ol):
        return ",".join([f"{o.oType}:{o.oId}" for o in ol])

    OG = ObjectGraph()
    a = OG.add("A", "o")
    b = OG.add("B", "1", deps=OG['A'])
    b = OG.add("B", "2", deps=OG['A'])
    b = OG.add("B", "3", deps=OG['A'])
    c = OG.add("C", "o", deps=OG['B'])

    assert s(a.deepDeps()) == ""
    assert s(b.deepDeps()) == "A:o"
    assert s(c.deepDeps()) == "B:1,B:2,B:3"
    assert s(c.deepDeps(level=1, mode='equal')) == "B:1,B:2,B:3"
    assert s(c.deepDeps(level=1, mode='lessOrEqual')) == "B:1,B:2,B:3"

    assert s(c.deepDeps(level=2)) == "A:o"
    assert s(c.deepDeps(level=2, mode='equal')) == "A:o"
    assert s(c.deepDeps(level=2, mode='lessOrEqual')) == "B:1,B:2,B:3,A:o"
    assert s(c.deepDeps(dot="A", level=2, mode='lessOrEqual')) == "A:o"
    assert s(c.deepDeps(dot="B", level=2, mode='lessOrEqual')) == "B:1,B:2,B:3"

    assert s(c.deepDeps(level=3)) == ""
    assert s(c.deepDeps(level=3, mode='equal')) == ""
    assert s(c.deepDeps(level=3, mode='lessOrEqual')) == "B:1,B:2,B:3,A:o"
    assert s(c.deepDeps(dot="A", level=3, mode='lessOrEqual')) == "A:o"
    assert s(c.deepDeps(dot="B", level=3, mode='lessOrEqual')) == "B:1,B:2,B:3"
