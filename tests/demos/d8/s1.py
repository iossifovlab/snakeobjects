#!/usr/bin/env python

from snakeobjects import Project

p = Project()


print('hi')
print('the project directory:', p.directory)
print('the pipeline directory:', p.get_pipeline_directory())
print('the parameters are:')
for k,v in p.parameters.items():
    print(f"\t{k}: {v}")
