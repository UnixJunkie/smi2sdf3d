#!/usr/bin/python

import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

if len (sys.argv) != 4:
    print "usage: %s N input.smi output.sdf" % sys.argv[0]
    sys.exit(1)

n_confs = int(sys.argv[1])
input_smi = sys.argv[2]
output_sdf = sys.argv[3]

reader = Chem.SmilesMolSupplier(input_smi)
writer = Chem.SDWriter(output_sdf)

def how_many_conformers(mol):
    nb_rot_bonds = AllChem.CalcNumRotatableBonds(mol)
    if nb_rot_bonds <= 7:
        return 50
    elif nb_rot_bonds <= 12:
        return 200
    else:
        return 300

rmsd_threshold = 0.35 # Angstrom

# keep only conformers which are far enough from the reference conformer
# (the one of lowest energy)
def rmsd_filter(mol, ref_conf, l):
    res = []
    refConfId = ref_conf.GetId()
    for e, curr_conf in l:
        currConfId = curr_conf.GetId()
        if AllChem.GetBestRMS(mol, mol, refConfId, currConfIf) > rmsd_threshold:
            res.append((e, curr_conf))
    return res

for mol in reader:
    if mol:
        n = how_many_conformers(mol)
        mol_H = Chem.addHs(mol)
        # FBR: check doc for EmbedMultipleConfs
        confIds = AllChem.EmbedMultipleConfs(mol_H, n)
        conf_energies = []
        # FF minimization
        for conf in confIds:
            confId = conf.get_Id()
            ff = AllChem.UFFGetMoleculeForceField(mol_h, conf_id=confId)
            ff.Minimize()
            energy = ff.CalcEnergy()
            conformer = mol_H.GetConformer(confId)
            conf_energies.append((energy, conformer))
        # remove neighbors
        # sort by increasing E
        sorted(conf_energies, key=lambda x: x[0])
        # output non neighbor ones
        nb_out = 0
        kept = []
        while nb_out < n_confs and len(conf_energies) > 0:
            ++nb_out
            (e, conf) = conf_energies.pop(0)
            keep.append((e, conf))
            # remove neighbors
            conf_energies = rmsd_filter(mol_h, conf, conf_energies)
        # write them out
        res = Chem.Mol(mol_h)
        res.RemoveAllConformers()
        for e, conf in kept:
            confId = conf.GetId()
            conf = mol_h.GetConformer(confId)
            res.AddConformer(conf, assignId=True)
        writer.write(res)
writer.close()

# FBR: check molecule names in output
