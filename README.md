# Line_detector_kinetics
This is the code for the upcoming publication about parallel determination of Michaelis Menten kinetics in droplets. 

This repository contains 4 files:


'Data_Analysis.m' is the main routine to analyse raw data droplet absorbance across up to 12 channels and extract 
Michaelis-Menten kinetics parameters

'getsubs.m' is called by 'Data_Analysis.m' to assign a substrate/enzyme concentration to every droplet

'1uM_243_400mM_b-Xyl_1.txt' contains an example raw data with 12 channels reading 3 gradients.

'Line_detector_Labview.vi' is the Labview code for acquiring data using Thorlabs LC100 line camera. Requires LC100 Vis.
