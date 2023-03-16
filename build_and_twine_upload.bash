#!/bin/bash

rm -rf package/dist/

python3 -m build package/

twine upload package/dist/*

