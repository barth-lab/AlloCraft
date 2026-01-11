# AlloCraft
AlloCraft is a two-step process that aims to modify allosteric behavior in proteins via amino acid design.

The first step of AlloCraft consists in predicting allosteric properties of your protein of interest by analyizing molecular dynamics (MD) simulations. Allosteric coupling is divided into a dynamic and structural element: dynamic allostery is infered by calculating mutual information (MI) for each pair of residues, followed by construction of residue networks (Pathways). The construction maximizes the MI along the path while minimizing its distance (i.e. number of residues connecting the distant pair). Constructed pathways are then clustered into allosteric pipelines. To infer allosteric structural perturbations, dihedral differences are calculated using Kullback-Liebler divergences between a reference and a target simulation. If the reference and target simulation correspond to an apo and ligand-bound states, respectively, then KLdiv will report the structural perturbations induced by ligand binding across the entire receptor. This allosteric analysis has been applied to GPCR-ligand systems, but could be applied to any other protein as well. The second step is the design step: design hotspots inferred from step 1 are mutated to all 20 possible amino acids using RosettaMembrane and assessing mutation stability using Rosetta scoring and coarse grained dynamical coupling using dynamic cross-correlation maps calculated using normal mode analysis. Final amino acid designs are chosen based on both stability and correlated dynamics scores.


![Slide1](https://github.com/barth-lab/AlloCraft/assets/68428149/37d137c3-b31a-415a-b990-fb7c0b2042a8)

## Step 1: Analysis of molecular dynamics simulations to extract allosteric properties:
This step offers 3 functionalities: **md2path**, **kldiv**, and **meta** analyses. Architecture of the code can be found in [reference.md](/step1_MD_simulation_analysis/reference.md)
- **md2path** is the main function that extracts allosteric pathways from MD trajectories of a given protein state using mutual information calculated from internal coordinates, all while providing a wide array of analysis methods along the way:
  1. Transform xtc to dcd (if needed)
  2. Fetch reference PDBs from database, and then align sequences and input PDBs to reference PDBs
  3. Load and align trajectories
  4. Calculate RMSD and RMSF of protein and ligand (if present)
  5. Calculate contact map with ligand and/or effector protein (if they exist)
  6. If effector protein is not present in the simulation but is present in reference PDB, contact map is calculated from reference PDB
  7. Calculate GPCR order parameters for activation states (if protein is a GPCR)
  8. PCA of ligand binding poses followed by clustering
  9. Calculate dihedral time series from trajectories (phi, psi, and chiX) 
  10. Calculate 1st and 2nd order entropies from dihedrals, followed by mutual information (MI)
  11. Assess convergence of entropies
  12. Run allosteric pathway calculation
  
- **kldiv** compares a test ensemble to a reference ensemble at a degree of freedom basis by calculating Kullback−Leibler divergences (KLdiv) between dihedrals:
  1. Fetch reference PDBs, and then align sequences and input PDBs to reference PDBs
  2. Load and align trajectories of test and reference ensembles
  3. Calculate dihedrals from trajectoriess (assumes same receptor IDs between entries)
  4. Reconcile dihedrals between test and reference systems (in case there are mutated residues with different dihedrals)
  5. Calculate 1st order KLdiv
  6. Calculate 2nd order KLdiv and mutual divergence (optional, and very slow, and not very informative)
  7. Plot dihedral distributions with highest divergences
  8. Visualize KLdivs residue-wise
  

- **meta** functions analyze the output of **md2path** or **kldiv** for a set of systems/proteins. Refer to individual scripts for more details

### Generated Data
The analysis pipeline generates a variety of data, which is typically saved in a md2pathdev subdirectory within your specific data folder (e.g., data/YOUR_SIMULATION_FOLDER/md2pathdev/). Key outputs include:

Plots and Figures (.fig, .pdf): Visualizations of RMSD, RMSF, contact maps, PCA results, and KL divergence.
PDB Files (.pdb): Structures with RMSF or KL divergence values saved in the B-factor column for easy visualization.
MAT-files (.mat): Contain processed data such as dihedral time series, MI matrices, and PCA results.
Text and Excel Files (.txt, .xlsx): Lists of interacting residues and other tabular data.

### Installation:

Requirements:
- [Matlab](https://www.mathworks.com/products/matlab.html) R2024a or newer
- [VMD](https://www.ks.uiuc.edu/Research/vmd/) [OPTIONAL]

Installation time: ~15 min (including the matlab installation)

The software has been tested with:
- Windows 11 version 21H2
- Ubuntu 18.04
- Matlab R2024a, R2024b, and R2025b [with a one-line fix]
- VMD 1.9.3

Download `step1_MD_simulation_analysis/`  folder either manually from this repo or by cloning it

Add the downloaded scripts and folders to the Matlab path
Adding to path in Matlab can be done in more than one way: 
- Per session:
https://ch.mathworks.com/help/matlab/ref/addpath.html
- Permanently at startup:
https://ch.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html

For the smoothest experience, AlloCraft also uses visual molecular dynamics (VMD) to transform .xtc files to .dcd and to predict secondary structure of a protein. You can download it for free from here: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD 

Don't forget to add VMD to your path variable too!

You can run the main scripts without having VMD installed if you provide .dcd trajectories instead of .xtc and if you provide *helices* variable in the input scripts (*input_md2path.m* and *input_kldiv.m*).

### Instructions:

Before running, you must configure two key files in the `step1_MD_simulation_analysis/` directory: `input_md2path.m` and `input_kldiv.m`. These scripts use a `global settings` variable to manage all parameters.

#### Important Parameters to Modify

##### 1. `input_md2path.m` (Main Analysis)

*   `settings.mydir`: **Crucial.** Set this to the directory containing your simulation data. The script expects this directory to contain subfolders like `run1`, `run2`, etc.
*   `settings.xtcName`: The name of your trajectory file within each `run*` folder. The default is `'traj.dcd'`. **Note:** The scripts are optimized for `.dcd` files. If you use `.xtc` files, a conversion will be attempted using a VMD Tcl script (`load_save.tcl`), but direct `.dcd` usage is recommended.
*   `settings.chains`: Defines the chain identifiers for the `[receptor, G protein, ligand]`. Use a hyphen `-` for any missing chains (e.g., `'A-B'` for a receptor and ligand). The default is `'ACB'`.
*   `settings.numRuns`: The total number of simulation runs (e.g., `run1`, `run2`, ...) to analyze from your data directory.
*   `settings.stride`: The step size for reading trajectory frames. A value of `1` reads every frame, while a value of `10` would read every 10th frame. This is useful for reducing computational load on very large trajectories.

##### 2. `input_kldiv.m` (Comparative Analysis)

*   `settings.mydir`: The data directory for your primary or "perturbed" system.
*   `settings.refdir`: **Reference System.** The data directory for the second simulation that you want to compare against the primary one. KL divergence measures the difference between the dihedral distributions of the `mydir` and `refdir` systems (e.g., comparing a mutant against a wild-type).
*   `settings.mainName` & `settings.refName`: Short, descriptive names for your main and reference systems used for labeling outputs.

##### 3. PDB Filenames

The scripts expect specific PDB filenames within the main analysis functions. If your files are named differently, you **must** update them.

*   In `md2pathMain.m`:
    ```matlab
    % Look for this line and change "prot.pdb" if needed
    database.read(fullfile(settings.mydir, "prot.pdb"), settings.chains, "Protein");
    ```
*   In `kldivMain.m`:
    ```matlab
    % Change "prot.pdb" in the reference or target directories if your files are named differently
    database.read(fullfile(settings.mydir, "prot.pdb"), settings.chains, settings.mainName);
    database.read(fullfile(settings.refdir, "prot.pdb"), settings.chains, settings.refName );
    ```

Create a PDB file with separate chains for receptor, ligand, and G-protein/effector site, and define the chains in the input script.

> **Note:** all input residues should be in the same pdb numbering as the input PDB.
#### Directory Structure:
Main directory:
prot.pdb (your pdb file)
- run1/ -> traj.dcd OR traj.xtc (your trajectory files)
- run2/ -> traj.dcd OR traj.xtc
- .
- .
- .
- runN/ -> traj.dcd OR traj.xtc

Where run1->N are directories that contain your trajectories

Run *input_md2path.m* or *input_kldiv.m*

Run *md2pathMain.m* or *kldivMain.m*. All the output can be found in a folder called **md2pathdev** created in your defined directory in your input script.

#### Typical runtime for calculations: 
- md2path: ~ 1 hour from start to finish on a typical PC and for around 15000-20000 frames of trajectory of a 300-400 residue GPCR
- kldiv: a few minutes assuming you only calculate 1st order KL divergences

### Running the Analysis

The analysis can be executed from the interactive MATLAB GUI.

1.  **Open MATLAB.**
2.  **Navigate to the scripts directory.** In the "Current Folder" toolbar, paste the absolute path to the `scripts2` directory and press Enter.
    ```matlab
    /path/to/AlloDy_Analysis_Repo/scripts2
    ```
3.  **Run the scripts.** In the "Command Window" (at the `>>` prompt), type the following commands. The semicolon at the end of each command suppresses verbose output and allows them to run sequentially.
    ```matlab
    input_md2path;
    md2pathMain;
    ```
    
    ```matlab
    input_kldiv;
    kldivMain;
    ```
    
    
### Acknowledgements:
- Mutual information statistical filtering: McClendon et al., J Chem Theory Comput (2009)
- Using Kullback−Leibler divergences to compare conformational ensembles: McClendon et al., J Chem Theory Comput (2012)
- Clustering mutual information into allosteric pathways: Bhattacharya and Vaidehi, Biophy J (2014); Nivedha et al., Mol Pharmacol (2018)
- MDToolbox for Matlab: Matsunaga, and Sugita, J Chem Phys (2018)

### FAQ:

#### 1- Compatibility with Matlab versions

| MATLAB | Status | Notes |
|--------|--------|-------|
| R2025b | ✅ OK | Needs 1-line fix (below). |
| R2024b | ✅ OK | No changes needed in our tests. |
| R2024a | ✅ OK | No changes needed in our tests. |


**Please include your MATLAB version when reporting issues.**

---

#### 2- Quick fix for R2025b

`calcPlotOrderParameters.m` may use `size(plots)` where `length` is safer.

**Path:** `scripts2/src/calcPlotOrderParameters.m`

**Change:**

```matlab
% Old
nplots = size(plots);

% New
nplots = length(plots);
```

**Then run:**

```matlab
input_md2path; md2pathMain; input_kldiv; kldivMain;
```

---

#### 3- R2024a first-run workaround

If `input_md2path.m` errors on the first run, try again.

If it persists:

```matlab
restoredefaultpath; rehash toolboxcache;
input_md2path; md2pathMain; input_kldiv; kldivMain;
```

## Step 2: Amino acid design to modify allosteric behavior
This is the second step of the process, where we will do in silico mutagenesis to the hubs that we extracted from step1, followed by coarse grained dynamical coupling calculations.

### Installation:
Requirements:
Rosetta software suit (https://www.rosettacommons.org/software/license-and-download)
Python (tested with python 3.6.7 and pandas 0.23.4)
R (bio3D package: http://thegrantlab.org/bio3d/)

### Instructions:
The first part is performing in silico mutagenesis to specific residues, which will be completed using RosettaMembrane:
For this step, we will need an input PDB and a span file that specifies the position of the membrane, for reference on how to prepare input files for membrane simulations with Rosetta:
https://www.rosettacommons.org/docs/latest/application_documentation/membrane_proteins/RosettaMP-GettingStarted-PreparingInputs#span-files

The python scripts provided will take care of running RosettaMembrane in an automated fashion.

We then extract the lowest energy pdbs from each mutagenesis and use as input for stability and dynamical coupling calculations.
The R scripts will perform NMA to extract DCCM and will extract energy differences upon mutation from the mutant PDB files.

More details can be found in the demo.

### Acknowledgements:
- bio3D from the Grant lab


## Contact: mahdi.hijazi@epfl.ch (2024)
##


