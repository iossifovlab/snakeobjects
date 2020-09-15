for p in proj*; do
    (cd $p; rm -rf *.snakefile OG.json OG.OG glbl.makefile log objLinks base)
done
