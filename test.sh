#!/bin/bash

\rm -f data/test_output.sdf

./smi2sdf.py -n 10 -i data/test_input.smi -o data/test_output.sdf
