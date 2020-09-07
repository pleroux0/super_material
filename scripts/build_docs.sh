#!/bin/sh

SPHINXBUILD=sphinx-build
SOURCEDIR="./docs/"
BUILDDIR="./dist/docs/"

#poetry run "${SPHINXBUILD}" "${SOURCEDIR} ""${BUILDDIR}"
poetry run "${SPHINXBUILD}" "${SOURCEDIR}" "${BUILDDIR}"
