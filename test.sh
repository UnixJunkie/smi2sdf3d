#!/bin/bash

out=data/test_output.sdf

\rm -f ${out}

./smi2sdf.py -n 10 -i data/test_input.smi -o ${out}

echo '#molecules in output:'
which molcount 2>&1 >/dev/null && molcount ${out}
