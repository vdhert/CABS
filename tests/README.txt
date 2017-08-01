To run test of choice type:

PYTHONPATH=`pwd` python tests/<file>.py -v

in cabsDock directory.

Some tests requires additional data stored in tests/data directory. If one would like to use another directory, use:

PYTHONPATH=`pwd` python tests/<file>.py -v --data_path my/data/path/
