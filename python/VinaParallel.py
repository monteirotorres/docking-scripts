#!/usr/bin/env python

import multiprocessing as mp
import subprocess
import os


def dock(ligand):
    try:
        ligandname = os.path.basename(ligand).replace(".pdbqt", "")
        outdir = os.path.join(os.getcwd(), ligandname)
        os.mkdir(outdir)
        outpath = os.path.join(outdir, "out.pdbqt")
        logpath = os.path.join(outdir, "log.txt")
        command = ["vina", "--config", "conf.txt", "--ligand", ligand, "--out", outpath]
        with open(logpath,'w') as f_obj:
            subprocess.run(command,stdout=f_obj,text=True)
        print(f"Done docking ligand {ligandname}; check output at: {outpath}")

    except subprocess.CalledProcessError:
        print(clrs['r']+'Command '+' '.join(command)+' returned non-zero status\n')
        pass


liganddir = input("Please state the path to the pdbqt ligands: ")

filelist = []
for file in os.listdir(liganddir):
    if file.endswith(".pdbqt"):
        filelist.append(os.path.join(liganddir, file))



p = mp.Pool(mp.cpu_count())

for result in p.imap_unordered(dock, filelist, chunksize=10):
    pass
