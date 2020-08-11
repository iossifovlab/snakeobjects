for p in proj*; do
    (cd $p; rm -rf .snakemake OG.json glbl.makefile log objLinks base)
done
