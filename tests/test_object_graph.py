import pytest
from iippl.ObjectGraph import ObjectGraph,ObjectGraphException

@pytest.fixture
def OG():
    return ObjectGraph()

def test_create(OG):
    assert OG.get_object_types() == [] 

def test_add(OG):
    o = OG.add('t','o')
    assert OG.get_object_types() == ['t'] 
    assert OG['t','o'] == o
    assert OG['t','o'].type == 't'
    assert OG['t','o'].name == 'o'
    assert OG['t','o'].deps == []
    assert OG['t','o'].params == {}
    assert OG['t'] == [o]

def test_dep_par(OG):
    OG.add('A','2')
    OG.add('A','1')
    OG.add('A','3')
    OG.add('B','o',deps=OG['A'])
    assert OG['B','o'].deps == [OG['A','2'],OG['A','1'],OG['A','3']]

def test_dep_par(OG):
    OG.add('t','o',params={'p':'v'})
    assert OG['t','o'].params == {'p':'v'}

def test_dep_par(OG):
    OG.add('t','o')
    with pytest.raises(ObjectGraphException) as excinfo:
        OG.add('t','o')
    assert "is already added" in str(excinfo.value)

def test_otype_order(OG):
    from random import random
    inOrder=[str(random()) for _ in range(100)]
    for t in inOrder:
        OG.add(t,'o')
    assert OG.get_object_types() == inOrder    
        

