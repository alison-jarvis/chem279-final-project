# CNDO/2 Molecular Property Classification

This repository contains code which runs a CNDO/2 based molecular classification pipeline on a set of structurally optimized input molecules. 

## Code Usage

### Overview

The code modules are contained in the `src` folder. `pipeline.cpp` handles running the actual pipeline, and all utilities corresponding to the CNDO/2 code and the dipole calculation are contained in `cndo2_utils.cpp`. Code pertaining to the approximated dipole calculation, polarizability, and derived properties is contained in `polarity.cpp`. 

### Container

This project uses container based development, with the already established `berkeley-chem-179-279/dev:latest` image. To launch into this container, run `interactive.sh`. 

### Input Data

The data corresponding to the molecules of interest are contained in the `input` directory. To add a molecule, you must create a .xyz file in the `atoms` directory, which contains the optimized structure of the molecule in question. You must also create a corresponding .json file in the `input` directory which references the name of the .xyz file, and the number of $\alpha$ and $\beta$ electrons. Additionally, you must ensure that the basis functions required for all atoms in a given molecule are contained within the `basis` folder. 

### Pipeline Run Instructions

Three main scripts allow you to run the pipeline. The script `build.sh` will compile the code. The script `run_pipeline.sh` will run the property classification pipeline for every molecule within the `input` folder. Finally, use the `clean.sh` script to remove previous outputs before re-building. All outputs from the pipeline will be written to a folder called `pipeline_outputs`. 

### Field Magnitude Run Instructions

Building will also create an executable called `test`, which corresponds to the test in which we varied the magnitude of the electric field and re-computed the polarizability results, to see if they were consistent across different field magnitude values. To perform this test, run `test.sh`. The results will be written out in csv file format to a folder `testing_outputs`. 
