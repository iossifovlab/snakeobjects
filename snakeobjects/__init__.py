"""
This is snakeobject python package.
"""
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .Project import Project
from .ObjectGraph import ObjectGraph, load_object_graph

