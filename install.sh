#! /bin/bash

# temporary, moze warto cos pozadnego napisac, ale moze nie

where='~/.local'
mode='develop'

python setup.py $mode --prefix=$where
