#!/usr/bin/env bash

sphinx-apidoc -fo docs/source bioflow
sphinx-build -b html docs/source docs/build >> doc_build.log 2>>doc_build_err.log
pandoc -o README.md inlined_readme.rst
