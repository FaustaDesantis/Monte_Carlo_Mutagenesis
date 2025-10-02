#!/bin/bash

# identification of the interfaces choosing 8Å as threshold for binding site definition
echo "************************INTERFACE IDENTIFICATION****************************"
Rscript Scripts/Identification_interfaces.R Data/xxxx_A.pdb Data/xxxx_B.pdb 8

# calculation of the initial descriptors considering chain A as fixed and chain B as mutable
echo "***********************INITIAL DESCRIPTORS CALCULATION**********************"
Rscript Scripts/Initial_descriptor_calculation.R Data/xxxx_A.pdb Data/xxxx_B.pdb Files/interface_prot_1.csv Files/interface_prot_2.csv

# Optimization protocol taking in input the initial descriptors
echo "***********************OPTIMIZATION PROTOCOL*************************"
Rscript Scripts/Monte_Carlo_optimization.R  Data/xxxx_B.pdb Files/interface_prot_2.csv Files/Surf_BS_1/xxxx_A_bs.csv Files/ZernDescr/Zernike_invariants_xxxx_A_bs.dat Files/PQR/xxxx_A.pqr Files/Initial_descriptors.csv
