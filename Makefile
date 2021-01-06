# Makefile for ClusterTrellis
SHELL := /bin/bash

# You can set these variables from the commandline.
VERSION=$(shell python setup.py --version)

./dist/ClusterTrellis-${VERSION}.whl:
	python setup.py sdist bdist_wheel

clean:
	rm dist/ClusterTrellis-${VERSION}.whl


install: ./dist/ClusterTrellis-${VERSION}-py3-none-any.whl # pip install
	pip install --upgrade ./dist/ClusterTrellis-${VERSION}.whl
	


%: Makefile