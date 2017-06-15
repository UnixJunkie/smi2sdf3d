#!/usr/bin/env python
# https://stackoverflow.com/questions/2429511/why-do-people-write-usr-bin-env-python-on-the-first-line-of-a-python-script

# generation of up to N low energy conformers
# from 2D input (smi) to 3D output (sdf)
# see Ebejer et. al.
# "Freely Available Conformer Generation Methods: How Good Are They?"
# JCIM, 2012, DOI: 10.1021/ci2004658 for technical details

# Copyright (C) 2017 Francois BERENGER
# System Cohort Division,
# Medical Institute of Bioregulation,
# Kyushu University
# 3-1-1 Maidashi, Higashi-ku, Fukuoka 812-8582, Japan

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Your code is Python 2 only because of the print statements. Use the import to replace it with the 
# print function and you have support for both (or even drop Python 2 completely)
from __future__ import print_function 

import sys
from contextlib import closing

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


# This is a Python 3 version of your code. It uses extended unpacking
def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            smile, name, *ignore = line.split()
            yield name, Chem.MolFromSmiles(smile)


# nb. conformers to generate prior to energy minimization
# as an empirical function of the molecule's flexibility
def how_many_conformers(mol):
    nb_rot_bonds = AllChem.CalcNumRotatableBonds(mol)
    if nb_rot_bonds <= 7:
        return 50
    elif nb_rot_bonds <= 12:
        return 200
    return 300  # This is more 

# keep only conformers which are far enough from the reference conformer
# (the one of lowest energy)
def rmsd_filter(mol, ref_conf, conf_energies):  # why not use the same name as in the calling code
    # I think the following is identical; of course you loose the ability for debug output
    refConfId = ref_conf.GetId()
    return [(e, curr_conf) for e, curr_conf in conf_energies
            if AllChem.GetBestRMS(mol, mol, refConfId, curr_conf.GetId()) > rmsd_threshold]
    # If I remember correctly, the GetBestRMS code is pretty slow if there are a lot of CH3 groups. 
    # You can get a considerable speedup by ignoring hydrogens. However, that is less trivial
#     # print("before: %d" % (len(l)))
#     res = []
#     for e, curr_conf in l:
#         currConfId = curr_conf.GetId()
#         rms = AllChem.GetBestRMS(mol, mol, refConfId, currConfId)
#         # print("e: %f rms: %f" % (e, rms))
#         if rms > rmsd_threshold:
#             res.append((e, curr_conf))
#     # print("after: %d" % (len(res)))
#     return res

# Separating the function definitions from the actual code is more a preference, but also makes it easier to 
# read the code

# I always use argparse for command line parameter handling. Gives a consistent interface. Probably won't even
# use more lines of code.
if len(sys.argv) != 4:
    print("usage: %s N input.smi output.sdf" % sys.argv[0])
    sys.exit(1)

n_confs = int(sys.argv[1])
input_smi = sys.argv[2]
output_sdf = sys.argv[3]

# to prune too similar conformers
rmsd_threshold = 0.35 # Angstrom

with closing(Chem.SDWriter(output_sdf)) as writer:
    for name, mol in RobustSmilesMolSupplier(input_smi):
        if mol is None:
            continue
# vv I would move this part into a separate function. 
        n = how_many_conformers(mol)
        print("init pool size for %s: %d" % (name, n))
        mol_H = Chem.AddHs(mol)
        print("generating starting conformers ...")
        conf_ids = AllChem.EmbedMultipleConfs(mol_H, n)  # Keep the naming style consistent
        # If the time for the embedding code is negligible, you can use the following for more concise code
        # for cid in AllChem.EmbedMultipleConfs(mol_H, n):
        conf_energies = []
        # FF minimization
        print("FF minimization ...")
        for cid in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(mol_H, confId=cid)
            # print("E before: %f" % ff.CalcEnergy())
            ff.Minimize()
            energy = ff.CalcEnergy()
            # print("E after: %f" % energy)
            conformer = mol_H.GetConformer(cid)
            # print("cid: %d e: %f" % (cid, energy))
            conf_energies.append((energy, conformer))
        # sort by increasing E
        conf_energies = sorted(conf_energies, key=lambda x: x[0])
        # output non neighbor conformers
        # nb_out = 0 # is always identical to len(kept)
        kept = []
        print("RMSD pruning ...")
        while len(kept) < n_confs and len(conf_energies) > 0:
            (e, conf) = conf_energies.pop(0)
            kept.append((e, conf))
            # remove neighbors 
            conf_energies = rmsd_filter(mol_H, conf, conf_energies)
# ^^ up to here
        # write them out
        print("kept %d confs for %s" % (len(kept), name))
        res = Chem.Mol(mol_H)
        res.RemoveAllConformers()
        for e, conf in kept:
            cid = res.AddConformer(conf, assignId=True)
            name_cid = "%s_%04d" % (name, cid)
            res.SetProp("_Name", name_cid)
            # print("cid: %d" % cid)
            writer.write(res, confId=cid)
