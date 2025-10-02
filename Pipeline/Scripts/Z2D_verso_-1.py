#!/usr/bin/env python
# coding: utf-8


import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl

import pandas as pd
from scipy.spatial import distance_matrix


sys.path.append("bin/")

import ZernikeFunc as ZF
import SurfaceFunc as SF

import time



res_path = "../Files/ZernDescr/"


Npixel = 25   # the plane edge in pixels...
Rs     = 6    # the radius of the sphere that includes the patch..
ZOrder = 20   # the Zernike expansion order..
VERSO = -1     # the orientation of the patch...


if(len(sys.argv) != 2):
	print("Error! Please insert:")
	print("python Z2D_verso_-1.py nome_file")
	print("nome_file, the file with the list of the surface")
	exit()

nome_file= sys.argv[1]
nome_folder =  "../Files/Surf_BS_2/"


## reading complex names and saving them in a list..
fp = open(nome_file)
mylist = []
for i in fp.readlines():
    mylist.append(i.split(".csv")[0])

    
for i in np.arange(0, len(mylist),1):   
        
        name = mylist[i]
        print(i)
        print("Processing protein %s\n"%(name))



      

        ## loading the surface patch obteined via the R script...
        
        mypatch_ = pd.read_csv("%s/%s.csv"%(nome_folder,name))

        
        l = len(mypatch_["x"])
        myl = len(mypatch_["x"])
        mypatch = np.zeros((l,6))

        mypatch[:,:] = mypatch_[["x","y","z", "Nx","Ny","Nz"]]
    
       
        surf = SF.Surface(mypatch[:,:3], patch_num = 0, r0 = Rs, theta_max = 45)

        
        rot_patch, rot_ag_patch_nv = surf.PatchReorientNew(mypatch, VERSO)
    

        # finding cone origin...                
        z = surf.FindOrigin(rot_patch)
        
           
        # creating plane.. 
        plane, weigths, dist_plane, thetas = surf.CreatePlane(patch=rot_patch, z_c=z , Np=Npixel)
        new_plane = surf.FillTheGap_everywhere(plane_=plane)

        ## enlarging plane..
        new_plane_ =  surf.EnlargePixels(new_plane)


        try:
            zernike_env.img  = new_plane_
        except:
            zernike_env = ZF.Zernike2d(new_plane_)

        coeff = zernike_env.ZernikeDecomposition(order=ZOrder)


        coeff_inv = np.absolute(coeff)


        np.savetxt("%s/Zernike_invariants_%s.dat"%(res_path, name), coeff_inv , fmt="%.6e")         

        



