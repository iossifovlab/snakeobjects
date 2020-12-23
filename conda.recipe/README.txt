# conda config --set anaconda_upload yes

# to add the gag
git tag -a 2.rc1 -m 'first 2 release candidate'
git push origin 2.rc1

Shold run: 
conda build . -c bioconda -c conda-forge


