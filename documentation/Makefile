# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

all: html 

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile html

helloWorld.tgz: 
	for d in `find helloWorld -name __pycache__`; do rm -r $$d; done
	for d in `find helloWorld -name objects`; do rm -r $$d; done
	tar czf $@ helloWorld

snakeobjectsTutorial.tgz: 
	for d in `find snakeobjectsTutorial -name __pycache__`; do rm -r $$d; done
	for d in `find snakeobjectsTutorial -name objects`; do rm -r $$d; done
	tar czf $@ snakeobjectsTutorial

snakeobjectsPaperExample.tgz: 
	for d in `find snakeobjectsPaperExample -name __pycache__`; do rm -r $$d; done
	for d in `find snakeobjectsPaperExample -name objects`; do rm -r $$d; done
	tar czf $@ snakeobjectsPaperExample

snakemakeTutorialExample.tgz: 
	for d in `find snakemakeTutorialExample -name __pycache__`; do rm -r $$d; done
	for d in `find snakemakeTutorialExample -name objects`; do rm -r $$d; done
	tar czf $@ snakemakeTutorialExample

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS). 
html: Makefile snakeobjectsTutorial.tgz snakeobjectsPaperExample.tgz snakemakeTutorialExample.tgz helloWorld.tgz
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

deploy: html
	rsync -avPHt -e "ssh -p 2022" _build/html lubo@iossifovlab.com:sites/snakeobjects/

clean:
	rm -rf _build
	rm -rf snakeobjectsTutorial.tgz
	rm -rf snakeobjectsPaperExample.tgz
	rm -rf snakemakeTutorialExample.tgz
