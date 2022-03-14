#!/usr/bin/env python

from snakeobjects import Project

p = Project()


print('hi')
print('the project directory:', p.directory)
print('the pipeline directory:', p.pipeline.get_definition())
print('the parameters are:')
for k, v in p.parameters.items():
    print(f"\t{k}: {v}")
