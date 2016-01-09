#!/usr/bin/env bash

sphinx-apidoc -fo docs/source src
sphinx-build -b html docs/source docs/build
