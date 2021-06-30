#!/usr/bin/env python
'''
<pipeline directory>
typeA.snakefile
Snakefile

<project directory> 				<----
	so_project.yaml
	.snakemake
		...
		...
	.sobjects
		OG.json

	objects			
		
		typeA
			O.A.1
				T.A.1.a     "typeA/O.A.1/T.A.1.a"
			O.A.2
		typeB
			O.B.1
'''
def get_remote_provider(provider): 
    if provider == "GS":
        from snakemake.remote.GS import RemoteProvider
    elif provider == "S3":
        from snakemake.remote.S3 import RemoteProvider
    else:
        raise Exception(f"Unknown provider {provider}")
    return RemoteProvider()

def get_project_files():
    # We may add additional project files through project configuratio option like, so_project_files=<glob 1>,<glob 2>,
    return [".sobjects/OG.json","so_projects.yaml"]

def upload_project_files_to_remote(provider,prefix):
    RP = get_remote_provider(provider)
    for f in get_project_files():
        RP.remote(f"{prefix}/{f}").upload_to_remote()

def download_project_files_from_remote(provider,prefix):
    RP = get_remote_provider(provider)
    for f in get_project_files():
        RP.remote(f"{prefix}/{f}").download_from_remote()


