#!/bin/sh

# Meep output to terminal
# NOTE: this file needs to be in 'MEEP code/src' to work

# Input parameters:
r1=.3      # melanosome radius
r2=0     # air radius / melanin radius
TE?=true   # cortex thickness
res=50      # resolution of simulation
nx=5        # number of melanosome layers
k2=.01      # extinction coeff for melanin

# Code to run meep
meep r=$r ra=$ra cor=$cor res=$res nx=$nx k2=$k2 TE?=true no-cortex?=false src/tri-rods-hollow_incid.ctl |tee struc0.out
meep r=$r ra=$ra cor=$cor res=$res nx=$nx k2=$k2 TE?=true no-cortex?=false no-struc?=false src/tri-rods-hollow_incid.ctl |tee struc1.out
grep flux1: struc0.out > flux0.dat
grep flux1: struc1.out > flux1.dat

# Plot output structure
h5topng output/eps-*.h5

# Plot reflectance spectrum using R
Rscript r_meep_plot.r
open Rplots.pdf
