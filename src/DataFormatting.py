# DataFormatting.py
# Python file for formatting a folder of PDB cif files into a pandas data 
# frame.

import os
import sys
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle

from rcsbsearchapi.search import TextQuery
from rcsbsearchapi import rcsb_attributes as attrs

def is_capital(string, space=''):
    for s in string:
        if s.islower() or (not s.isupper() and s!=space):
            return False
    return True

def atom_counts(f):
    with open(f, "r") as outfile:
        data = outfile.readlines()
    counts = {}
    for line in data:
        if 'ATOM' not in line:
            continue
        split = line.split()
        if split[2] in counts.keys():
            counts[split[2]] = counts[split[2]]+1
        elif split[2].isalpha() and len(split[2])<3 and split[2].isupper():
            counts.update({split[2] : 1})
    return counts

def acid_counts(f):
    with open(f, "r") as outfile:
        data = outfile.readlines()
    for i in range(len(data)):
        if "_entity_poly_seq.entity_id " in data[i]:
            n = i
            break
    m = data[n:].index("# \n")
    counts = {}
    for line in data[n+1:n+m]:
        split = line.split()
        if not (split[0].isnumeric() and split[1].isnumeric()):
            continue
        if split[2] in counts.keys():
            counts[split[2]] = counts[split[2]]+1
        elif len(split[2])==3 and split[2].isalpha() and split[2].isupper():
            counts.update({split[2] : 1})
    return counts

def bond_counts(f):
    with open(f, "r") as outfile:
        data = outfile.readlines()
    n = -1
    for i in range(len(data)):
        if "_chem_comp_bond.comp_id " in data[i]:
            n = i
            break
    if n == -1:
        return {}
    m = data[n:].index("# \n")
    counts = {}
    for line in data[n+1:n+m]:
        split = line.split()
        if not is_capital(split[0]):
            continue
        if split[3] in counts.keys():
            counts[split[3]] = counts[split[3]]+1
        elif split[3].islower():
            counts.update({split[3] : 1})
    return counts

def oligomeric_count(f):
    with open(f, "r") as outfile:
        data = outfile.readlines()
    n = -1
    for i in range(len(data)):
        if "_pdbx_struct_assembly.oligomeric_count" in data[i]:
            n = i
            break
    if len(data[i].split())>1:
        return int(data[i].split()[1])
    else:
        m = data[n:].index("# \n")
        for line in data[n:n+m]:
            split = line.split()
            if len(split)>1:
                return int(split[-1])

def protein_data_frame(source_dir, target=None):

    location = os.path.join(source_dir, '*.cif')
    filenames = glob.glob(location)

    all_atom_counts = {}
    all_acid_counts = {}
    all_bond_counts = {}
    index = []
    n = 0
    for f in filenames:
        n = n+1
        txt = (("\rProcessing protein {current} of {total}." 
                + "   {percent}% complete.")
                .format(current = n, total = len(filenames), 
                        percent = math.floor(n/len(filenames)*100)))
        sys.stdout.write(txt)
        sys.stdout.flush()
        
        all_atom_counts.update({f : atom_counts(f)})
        all_acid_counts.update({f : acid_counts(f)})
        all_bond_counts.update({f : bond_counts(f)})
        with open(f, "r") as outfile:
            data = outfile.readlines()
        for line in data:
            if "_entry.id" in line:
                index.append(line.split()[1])
                break
    all_atom_types = set().union(*[set(counts.keys()) for counts in all_atom_counts.values()])
    all_acid_types = set().union(*[set(counts.keys()) for counts in all_acid_counts.values()])
    all_bond_types = set().union(*[set(counts.keys()) for counts in all_bond_counts.values()])
    
    data = {}
    for atom_type in all_atom_types:
        data.update({atom_type : [all_atom_counts[f].get(atom_type, 0) for f in filenames]})
    for acid_type in all_acid_types:
        data.update({acid_type : [all_acid_counts[f].get(acid_type, 0) for f in filenames]})
    for bond_type in all_bond_types:
        data.update({bond_type : [all_bond_counts[f].get(bond_type, 0) for f in filenames]})
    data.update({"OLIG_COUNT" : [oligomeric_count(f) for f in filenames]})
    with open("symmetry-lists/C2_list.pkl", "rb") as file:
        symmetry_list = pickle.load(file)
    data.update({"SYMMETRY" : [int(f[-8:-4].upper()+'-1' in symmetry_list) for f in filenames]})

    if target is not None:
        data.to_csv(target)
    
    return pd.DataFrame(data, index=index)
        