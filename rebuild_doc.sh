#!/usr/bin/env bash

sphinx-apidoc -fo docs/source src
sphinx-build docs/source docs/build
