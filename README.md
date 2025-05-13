# HSES Evaluation Report

## Project Overview

This report provides a comprehensive evaluation of the HSES simulation setup, including compilation flags, default configuration, and plotting scripts. It aims to help users verify that the system remains in hydrostatic equilibrium.

---

## Compilation Flags

**Enabled**: `MODEL=HYDRO`  
**Disabled**: `GRAVITY`, `PARTICLE`

---

## Default Setup

1. Adopt Lohner’s error estimator on pressure as the refinement criterion  
   → Refinement threshold: `0.80`

2. The maximum refinement level (`MAX_LEVEL`) is set to `2` in the central region, while the outer region remains at level `0`.

---

## Notes

1. The following 4 plotting scripts are provided:

   - `plot_Press.py`  
   - `plot_Rho.py`  
   - `plot_slice_gas.py`  
     → Uses yt’s `SlicePlot` to create a density slice through the center of the simulation box.  
   - `radial_profile_evolution.py`  
     → Uses yt’s `ProfilePlot` to plot the time evolution of the radial density profile and compare it with the initial analytical density, thereby verifying that the system remains in hydrostatic equilibrium.
   - `mass_conservation.py`
     → Check mass status  

2. For the derivation and explanation of the pressure equation, please refer to the HackMD document:  
    [https://hackmd.io/mH2qiL4zRii5Pbz6Tn6ZcA?view](https://hackmd.io/mH2qiL4zRii5Pbz6Tn6ZcA?view)
