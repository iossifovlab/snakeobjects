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

snakeobjectsTutorial.tgz: 
	for d in `find snakeobjectsTutorial -name __pycache__`; do rm -r $$d; done
	for d in `find snakeobjectsTutorial -name objects`; do rm -r $$d; done
	tar czf $@ snakeobjectsTutorial

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
html: Makefile snakeobjectsTutorial.tgz 
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

deploy: html
	rsync -avPHt -e "ssh -p 2022" _build/html lubo@iossifovlab.com:sites/snakeobjects/

clean:
	rm -rf _build
	rm -rf snakeobjectsTutorial.tgz
