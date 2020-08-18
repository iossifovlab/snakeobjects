for p in proj*; do
    (cd $p; rm -rf .snakemake OG.json OG.OG glbl.makefile log objLinks base)
done
