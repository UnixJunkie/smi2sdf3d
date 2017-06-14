#!/usr/bin/python

import sys

if len (sys.argv) != 4:
    print "usage: %s N input.smi output.sdf" % sys.argv[0]
    sys.exit(1)

n_confs = int(sys.argv[1])
input_smi_fn = sys.argv[2]
output_sdf_fn = sys.argv[3]
