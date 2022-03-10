#!/usr/bin/env python3
# ----------------------------------------------------------
# Copyright (C) 2017 PHARAMACELERA S.L.
# All rights reserved.
# 
# WARRANTY DISCLAIMER
#
# THESE MATERIALS ARE PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PHARMACELERA OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THESE
# MATERIALS, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# File: genConf.py
#
# Created on 19/07/2017
# Last update 29/11/2021
# ----------------------------------------------------------
#
# Molecular conformer generator
#
# Example:
#
# genConf.py -i file_input.sdf -o file_output.sdf
# -n number_of_conformers (optional, if not specified is based
# on the nomber of rotable bonds) -rms rms_threshold
# -e energy_window (optional, Kcal/mol) -t number_of_threads (if not specify 1)
# -ETKDG (optional, use the useExpTorsionAnglePrefs and useBasicKnowledge options)
# -logfile (Default: False, write a log file with info about the conformer generation)
# ----------------------------------------------------------
# ----------------------------------------------------------
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from concurrent import futures
import progressbar
import argparse
import copy
import time
import gzip

def genConf(m, nc, rms, efilter, rmspost, nmol, molName, etkdg):
    mini=None
    info = ""
    nr = int(AllChem.CalcNumRotatableBonds(m))
    m4 = copy.deepcopy(m)
    m = Chem.AddHs(m)
    Chem.AssignAtomChiralTagsFromStructure(m, replaceExistingTags=True)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    doublestereolist=[]
    for b in range(0, m.GetNumBonds()):
        st = str(m.GetBondWithIdx(b).GetStereo())
        ty = str(m.GetBondWithIdx(b).GetBondType())
        if ty == "DOUBLE" and (st == 'STEREOZ' or st == "STEREOCIS"):
                doublestereolist.append((b, "STEREO Z"))
        elif ty == "DOUBLE" and (st == "STEREOE" or st == "STEREOTRANS"):
                doublestereolist.append((b, "STEREO E"))
    Chem.AssignStereochemistry(m, cleanIt=True, flagPossibleStereoCenters=True)
    chiralcenter = Chem.FindMolChiralCenters(m)+doublestereolist
    if nc == "X":
        if nr < 3:
            nc = 50
        elif nr > 6:
            nc = 300
        else:
            nc = nr**3
    m3 = copy.deepcopy(m)
    m5 = copy.deepcopy(m)
    ids=AllChem.EmbedMultipleConfs(m, numConfs=nc, pruneRmsThresh=rms, randomSeed=1, useExpTorsionAnglePrefs=etkdg, useBasicKnowledge=etkdg)
    if rms == -1 and efilter == "Y":
        if len(ids) != nc:
            info = "WARNING: " + molName + " generated less molecules than those required\n"
    numconf = m.GetNumConformers()
    if numconf == 0:
        m = copy.deepcopy(m3)
        ids=AllChem.EmbedMultipleConfs(m, numConfs=nc, pruneRmsThresh=rms, randomSeed=1)
        info = "WARNING: Molecule number " + str(nmol) + " embed without ETKDG method, molecule name: " + molName + "\n"
    m2 = copy.deepcopy(m)
    diz = []
    diz2 = []
    diz3 = []
    if m.GetNumConformers() == 0:
        info = "ERROR: Impossible to generate conformers of molecule " + str(nmol) + ", molecule name: " + molName + "\n"
        o = m4
        o = Chem.AddHs(o)
        embd = AllChem.EmbedMolecule(o, randomSeed=1)
        if embd == -1:
            info = "ERROR: Impossible to generate conformers of molecule " + str(nmol) + ", molecule name: " + molName + "\n"
            o = m4
        diz3 = [(None, -1)]

    else:
        for id in ids:
            molec = m.GetConformer(id).GetOwningMol()
            doublestereolist=[]
            for b in range(0, molec.GetNumBonds()):
                st = str(molec.GetBondWithIdx(b).GetStereo())
                ty = str(molec.GetBondWithIdx(b).GetBondType())
                if ty == "DOUBLE" and (st == 'STEREOZ' or st == "STEREOCIS"):
                    doublestereolist.append((b, "STEREO Z"))
                elif ty == "DOUBLE" and (st == "STEREOE" or st == "STEREOTRANS"):
                    doublestereolist.append((b, "STEREO E"))
            Chem.AssignStereochemistry(molec, cleanIt=True, flagPossibleStereoCenters=True)
            confchiralcenter = Chem.FindMolChiralCenters(molec)+doublestereolist
            if confchiralcenter != chiralcenter:
                m.RemoveConformer(id)
        if m.GetNumConformers() == 0:
            m = m5
            ids=AllChem.EmbedMultipleConfs(m, numConfs=nc, pruneRmsThresh=rms, randomSeed=1)
            for id in ids:
                molec = m.GetConformer(id).GetOwningMol()
                for b in range(0, molec.GetNumBonds()):
                    st = str(molec.GetBondWithIdx(b).GetStereo())
                    ty = str(molec.GetBondWithIdx(b).GetBondType())
                    if ty == "DOUBLE" and (st == 'STEREOZ' or st == "STEREOCIS"):
                        doublestereolist.append((b, "STEREO Z"))
                    elif ty == "DOUBLE" and (st == "STEREOE" or st == "STEREOTRANS"):
                        doublestereolist.append((b, "STEREO E"))
                Chem.AssignStereochemistry(molec, cleanIt=True, flagPossibleStereoCenters=True)
                confchiralcenter = Chem.FindMolChiralCenters(molec)+doublestereolist
                if confchiralcenter != chiralcenter:
                    if info != "":
                        info = info + "\n"
                    info = "WARNING: one or more conformer of Molecule number " + str(nmol) + " were excluded becouse it/they has/have different isomerism respect the input: " + molName + "\n"
                    m.RemoveConformer(id)
        try:
            if AllChem.MMFFHasAllMoleculeParams(m) == True:
                sm = copy.deepcopy(m)
                try:
                    for id in ids:
                        prop = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
                        ff = AllChem.MMFFGetMoleculeForceField(m, prop, confId=id)
                        ff.Minimize()
                        en = float(ff.CalcEnergy())
                        econf = (en, id)
                        diz.append(econf)
                except:
                    m = sm
                    for id in ids:
                        ff = AllChem.UFFGetMoleculeForceField(m, confId=id)
                        ff.Minimize()
                        en = float(ff.CalcEnergy())
                        econf = (en, id)
                        diz.append(econf)
                    if info != "":
                        info = info + "WARNING: Molecule number " + str(nmol) + " optimized with UFF force field, molecule name: " + molName + "\n"
                    else:
                        info = "WARNING: Molecule number " + str(nmol) + " optimized with UFF force field, molecule name: " + molName + "\n"
            else:
                for id in ids:
                    ff = AllChem.UFFGetMoleculeForceField(m, confId=id)
                    ff.Minimize()
                    en = float(ff.CalcEnergy())
                    econf = (en, id)
                    diz.append(econf)
                if info != "":
                    info = info + "WARNING: Molecule number " + str(nmol) + " optimized with UFF force field, molecule name: " + molName + "\n"
                else:
                    info = "WARNING: Molecule number " + str(nmol) + " optimized with UFF force field, molecule name: " + molName + "\n"
        except:
            m = m2
            if info != "":
                info = info + "ERROR: Unable to minimize molecule number: " + str(nmol)+ ", molecule name: " + molName + "\n"
            else:
                info = "ERROR: Unable to minimize molecule number: " + str(nmol)+ ", molecule name: " + molName + "\n"
            for id in ids:
                en = None
                econf = (en, id)
                diz.append(econf)

        if efilter != "Y":
            n, diz2, mini = ecleaning(m, diz, efilter)
        else:
            n = m
            diz2 = copy.deepcopy(diz)
            diz.sort()
            mini = float(diz[0][0])
            mini
        if rmspost != False and n.GetNumConformers() > 1:
            o, diz3 = postrmsd(n, diz2, rmspost)
        else:
            o = n
            diz3 = diz2

    return o, diz3, nr, info, mini

def ecleaning(m, diz, efilter):
    diz.sort()
    mini = float(diz[0][0])
    sup = mini + efilter
    n = Chem.Mol(m)
    n.RemoveAllConformers()
    n.AddConformer(m.GetConformer(int(diz[0][1])))
    diz2=[[float(diz[0][0]), int(diz[0][1])]]
    del diz[0]
    for x,y in diz:
        if x <= sup:
            n.AddConformer(m.GetConformer(int(y)))
            uni = [float(x), int(y)]
            diz2.append(uni)
        else:
            break
    return n, diz2, mini

def postrmsd(n, diz2, rmspost):
    diz2.sort()
    o = Chem.Mol(n)
    o.RemoveAllConformers()
    confidlist = [diz2[0][1]]
    enval = [diz2[0][0]]
    nh = Chem.RemoveHs(n)
    nh = Chem.DeleteSubstructs(nh, Chem.MolFromSmiles('F'))
    nh = Chem.DeleteSubstructs(nh, Chem.MolFromSmiles('Br'))
    nh = Chem.DeleteSubstructs(nh, Chem.MolFromSmiles('Cl'))
    nh = Chem.DeleteSubstructs(nh, Chem.MolFromSmiles('I'))
    del diz2[0]
    for z,w in diz2:
        confid = int(w)
        p=0
        for conf2id in confidlist:
            rmsd = AllChem.GetBestRMS(nh, nh, prbId=confid, refId=conf2id)
            if rmsd < rmspost:
                p=p+1
                break
        if p == 0:
            confidlist.append(int(confid))
            enval.append(float(z))
    for id in confidlist:
        o.AddConformer(n.GetConformer(id))
    diz3 = zip(enval, confidlist)
    return o, diz3

#----------------------------------- Main function --------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Molecular conformer generator')
    parser.add_argument('-i', required=True, help='sdf input file')
    parser.add_argument('-o', required=True, help='sdf output file')
    parser.add_argument('-n', type=int, required=False, help='number of conformers')
    parser.add_argument('-rms', type=float, required=False, help='rms threshold pre optimization')
    parser.add_argument('-e', type=float, required=False, help='energy window')
    parser.add_argument('-notprintproperty', action='store_true', help='Print molecule properties (energy and rotable bond number)')
    parser.add_argument('-ETKDG', action='store_true', help='use the useExpTorsionAnglePrefs and useBasicKnowledge options')
    parser.add_argument('-logfile', action='store_true', help='write a log file')
    parser.add_argument('-t', type=int, required=False, help='number of threads')
    args = parser.parse_args()
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    time_start = time.time()
    inp = args.i
    out = args.o

    if args.rms == None:
        rmspre = -1
    else:
        rmspre = args.rms
    if args.rms == None:
        rmspost = False
        irmsd = "Not"
    else:
        rmspost = args.rms
        irmsd = str(args.rms) + " A"
    if args.t == None:
        threads = 1
    else:
        threads = args.t
    if args.e == None:
        ent = "Y"
        iener = "Not"
    else:
        ent = args.e
        iener = str(args.e)+" Kcal/mol"
    if args.n==None:
        nc = "X"
        inc = "Rotatable bond based"
    else:
        nc = args.n
        inc = str(args.n)
    etkdg = False
    etkdginfo = "Not"
    if args.ETKDG == True:
        etkdg = True
        etkdginfo = "Yes"

    writer = Chem.SDWriter(out)

    if inp.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(inp)
    elif inp.endswith(".gz"):
        inp = gzip.open(inp)
        suppl = Chem.ForwardSDMolSupplier(inp)

    totmol = len(suppl)

##log file##
    if args.logfile == True:
        log = open(out.replace(".sdf", ".log"), "w")
        log.write("Conformer generator\nTotal molecule number: " + str(totmol) + "\nEnergy filter: " + iener + "\nRMSD filter: " + irmsd + "\nETKDG method: " + etkdginfo + "\nNumber of conformers per molecule: " + inc + "\nNumber of threads: "+str(threads)+"\n\n")
    ncicle = (totmol//1000)
    numMol = 0
    NM = []
    with futures.ProcessPoolExecutor(max_workers=threads) as executor:
        start = 0
        if totmol <= 999:
            end = totmol-1
        else:
            end = 999
        for cicle in range(0, ncicle+1):
            info = ""
            jobs = []
            nm = []
            for sup in range(start,end+1):
                numMol = numMol+1
                mol = suppl[sup]
                if mol is not None:
                    NM.append(numMol)
                    molName = mol.GetProp('_Name')
                    nm.append(molName)
                    job = executor.submit(genConf, mol, nc, rmspre, ent, rmspost, numMol, molName, etkdg)
                    jobs.append(job)
                else:
                    info = "ERROR: Impossible read " + str(numMol)
                    if args.logfile == True:
                        log.write(info+ "\n")
            widgets = ["Generating conformations of molecules from " + str(start+1) + " to " + str(end+1) +"; " , progressbar.Percentage(), " ", progressbar.ETA(), " ", progressbar.Bar()]
            pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(jobs))
            for job in pbar(futures.as_completed(jobs)):
                mol,ids,nr,info,mini = job.result()
            for j in range(0, len(jobs)):
                mol,ids,nr,info,mini = jobs[j].result()
                name = nm[j]
                propDic = suppl[NM[j+start]-1].GetPropsAsDict()
                if info != "":
                    mol.SetProp('genConf_Info', info)
                    if args.logfile == True:
                         log.write(info)
                for en,id in ids:
                    mol.SetProp('_Name', name)
                    if args.notprintproperty == False:
                        for k in propDic:
                            mol.SetProp(k, str(propDic[k]))
                        if en == None:
                            mol.SetProp('CONF_ENERGY', "Unable to calculate")
                        else:
                            mol.SetProp('CONF_ENERGY', str(en)+ " Kcal/mol")
                            if mini != en:
                                mol.SetProp('ENERGY_DIFFERENCE', str(en-mini)+ " Kcal/mol")
                            else:
                                mol.SetProp('ENERGY_DIFFERENCE', "0.00 Kcal/mol")
                        mol.SetProp('RotatableBondsNumber', str(nr))
                        if args.logfile == False and info != "":
                            print(info.replace('\n',''))
                        if not "ERROR: Impossible to generate conformers of molecule " in info:
                            writer.write(mol, confId=id)
            start = start + 1000
            if cicle == ncicle -1:
                end = totmol-1
            else:
                end = end + 1000
            if end == sup:
                break
    writer.close()
    time_end = time.time()
    exTime = time_end - time_start
    if args.logfile == True:
        log.write("Execution time " + str(exTime) + " seconds")
        log.close()
