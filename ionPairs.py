#!/usr/bin/env python2
# This script will take a specified .pdb file and estimate the number of ion-pairs.
# This number is based on side-chain oxygen atoms of Glu or Asp within a specified
# distance of any nitrogen atom from His, Arg or Lys.
# The .pdb file is downloaded into local dir "./pdbFiles".
# The script assumes modules Bio and Bio.PDB from Biopython are installed and in your python-path.
#
# Use at own discretion. It is very likely that a more elegant solution could be or have been written. 
#
# Author: Tobias Tandrup
# Date: 20200121
# Institution: UCPH

import Bio
import Bio.PDB as PDB
from Bio.PDB import PDBParser
from Bio.PDB import PDBList
pdb_parser = PDBParser(QUIET=True)


# User input PDB file and define ion-pair distance
pdb_file = raw_input("Enter a PDB entry:  ")
ion_pair_distance = float(raw_input("How far is an ion-pair distance [Angstrom]?  "))


def oxy_to_nit(res_oxy, res_nit, number_of_ion_pairs, ion_pair_distance):
# Define the relevant atoms and update the number of ion-pairs
    if res_nit.get_resname() == "HIS":
        # If the HIS residue is encountered in the sequence,
        # measure the distance between N atoms and the current oxygen atom.
        distance = res_oxy[oxy] - res_nit["ND1"]
        if distance <= ion_pair_distance:
            number_of_ion_pairs_changed = 1
            return number_of_ion_pairs_changed
        distance = res_oxy[oxy] - res_nit["NE2"]
        if distance <= ion_pair_distance:
            number_of_ion_pairs_changed = 1
            return number_of_ion_pairs_changed
    if res_nit.get_resname() == "ARG":
        # Same for ARG.
        distance = res_oxy[oxy] - res_nit["NH1"]
        if distance <= ion_pair_distance:
            number_of_ion_pairs_changed = 1
            return number_of_ion_pairs_changed
        distance = res_oxy[oxy] - res_nit["NH2"]
        if distance <= ion_pair_distance:
            number_of_ion_pairs_changed = 1
            return number_of_ion_pairs_changed
        distance = res_oxy[oxy] - res_nit["NE"]
        if distance <= ion_pair_distance:
            number_of_ion_pairs_changed = 1
            return number_of_ion_pairs_changed
    if res_nit.get_resname() == "LYS":
        # Same for Lys.
        distance = res_oxy[oxy] - res_nit["NZ"]
        if distance <= ion_pair_distance:
            number_of_ion_pairs_changed = 1
            return number_of_ion_pairs_changed
    return 0

# Define relevant N containing residues.
list_of_nit = ["HIS", "ARG", "LYS"]
# Define relevant O containing residues.
list_of_oxy = ["ASP", "GLU"]

# Import the PDB
PDBList().retrieve_pdb_file(pdb_file, pdir="./pdbFiles/")
name = pdb_parser.get_structure("name", "./pdbFiles/pdb" + pdb_file + ".ent")

# Go through the structure, find the relevant residues,
number_of_ion_pairs = 0
temp_number_of_ion_pairs = 0
res_oxy = []
res_nit = []
for res in name.get_residues():
    if res.get_resname() in list_of_oxy:
        res_oxy.append(res)
    if res.get_resname() in list_of_nit:
        res_nit.append(res)


# Check for ion pairs between O containing and N containing residues.
for i in range(len(res_oxy)):
    if res_oxy[i].get_resname() == "GLU":
        temp_number_of_ion_pairs = number_of_ion_pairs
        oxy = "OE1"
        for j in range(len(res_nit)):
            # Update number of ion-pairs.
            number_of_ion_pairs_changed = oxy_to_nit(res_oxy[i],
                                                     res_nit[j],
                                                     number_of_ion_pairs,
                                                     ion_pair_distance)

            number_of_ion_pairs = number_of_ion_pairs + number_of_ion_pairs_changed


        if temp_number_of_ion_pairs == number_of_ion_pairs:
            oxy = "OE2"
            for j in range(len(res_nit)):
                number_of_ion_pairs_changed = oxy_to_nit(res_oxy[i],
                                                         res_nit[j],
                                                         number_of_ion_pairs,
                                                         ion_pair_distance)

                number_of_ion_pairs = number_of_ion_pairs + number_of_ion_pairs_changed

    if res_oxy[i].get_resname() == "ASP":
        temp_number_of_ion_pairs = number_of_ion_pairs
        oxy = "OD1"
        for j in range(len(res_nit)):
            number_of_ion_pairs_changed = oxy_to_nit(res_oxy[i],
                                                     res_nit[j],
                                                     number_of_ion_pairs,
                                                     ion_pair_distance)

            number_of_ion_pairs = number_of_ion_pairs + number_of_ion_pairs_changed


        if temp_number_of_ion_pairs == number_of_ion_pairs:
            oxy = "OD2"
            for j in range(len(res_nit)):
                number_of_ion_pairs_changed = oxy_to_nit(res_oxy[i],
                                                         res_nit[j],
                                                         number_of_ion_pairs,
                                                         ion_pair_distance)

                number_of_ion_pairs = number_of_ion_pairs + number_of_ion_pairs_changed

    # Print output.
    print "\n\n>> Number of Ion-pairs in PDB={} is {}, based on a cut-off of {} Angstrom.\n".format(pdb_file, number_of_ion_pairs, ion_pair_distance)
