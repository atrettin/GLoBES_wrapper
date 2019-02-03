#!/usr/bin/env python 

import os, sys

import numpy as np


globes_path = "/afs/ifh.de/user/t/terliuk/scratch/software/GLoBES_test_inst/GLoBES_wrapper/wrapper/"
sys.path.append(globes_path)
import GLoBES

PREM_model = np.array( [
    [    0.00, 13.0],
    [1221.50, 13.0],
    [1974.33, 11.948163],
    [2727.17, 11.341572],
    [3480.00, 10.462371],
    [4220.33, 5.381883],
    [4960.67, 5.007959],
    [5701.00, 4.603493],
    [5971.00, 3.9],
    [6346.60, 3.4],
    [6356.00, 2.900000],
    [6368.00, 2.600000],
    [6371.00, 1.020000]
    ])

curdir = os.getcwd()
os.chdir(globes_path)
globes_calc =  GLoBES.GLoBESCalculator("test")
os.chdir(curdir)
globes_calc.InitSteriles(2)
globes_calc.SetEarthModel(PREM_model[:,0].tolist(), PREM_model[:,1].tolist())
globes_calc.PrintEarthModel()


DM21   = 7.60e-5
DM31   = 2.35e-3 + 0.5*DM21
DM41 = 1.0

sin2_23 = 0.5
sin2_13 = 0.0235
sin2_12 = 0.312
sin2_a = 0.04

theta23 = np.arcsin(np.sqrt(sin2_23))
theta13 = np.arcsin(np.sqrt(sin2_13))
theta12 = np.arcsin(np.sqrt(sin2_12))

deltacp = 0.

theta14 = 0.0
theta24 = np.arcsin(np.sqrt(0.02))
theta34 = np.arcsin(np.sqrt(0.0))

params = np.array( [theta12, theta13, theta23, 
                    deltacp, DM21, DM31, 
                    DM41, theta14,  theta24, theta34, 0.0, 0.0 ], dtype=float )

globes_calc.SetParametersArr(params ) 
globes_calc.PrintParameters()

energies = 10.0**np.linspace(0.0,4.0,301)
costheta = -1.0*np.ones_like(energies)
flavin   = -14*np.ones_like(energies)
flavout  = -14*np.ones_like(energies)

prob_numubar_numubar_sterile = globes_calc.MatterProbPDGArr(flavin,
                             flavout, 
                             energies, 
                             costheta)

params = np.array( [theta12, theta13, theta23, 
                    deltacp, DM21, DM31, 
                    DM41, theta14,  0, 0, 0.0, 0.0 ], dtype=float )
globes_calc.SetParametersArr(params ) 

prob_numubar_numubar_std = globes_calc.MatterProbPDGArr(flavin,
                             flavout, 
                             energies, 
                             costheta)
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))
plt.plot(energies, prob_numubar_numubar_sterile, label = "sterile" )
plt.plot(energies, prob_numubar_numubar_std, label = "standard")

plt.xscale("log")
plt.ylim(0,1)
plt.legend()
plt.xlabel("Energy [ GeV ]")
plt.ylabel("Probability")
plt.savefig("probability.pdf")
