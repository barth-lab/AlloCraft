# AlloDy
A package to analyze molecular dynamics trajectories of proteins, with a focus on allosteric pathways and ensemble differences. Description of the package, installation, and instructions for running can be found here. Architecture of the code can be found in [reference.md](/reference.md)


The package offers 3 functionalities: **md2path**, **kldiv**, and **meta** analyses. 
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

## Installation:

Requirements:
- [Matlab](https://www.mathworks.com/products/matlab.html) 2021b or newer
- [Matlab Bioinformatics toolbox](https://www.mathworks.com/products/bioinfo.html)
- [MDprot](https://github.com/mahdiofhijaz/MDprot) package
- [VMD](https://www.ks.uiuc.edu/Research/vmd/)

The software has been tested with:
- Windows 11 version 21H2
- Ubuntu 18.04
- Matlab R2021b
- VMD 1.9.3

Download MDprot and add it to your Matlab path.
Downloading can be done manually from the aforementioned link or in the command line via:

```
git clone https://github.com/mahdiofhijaz/MDprot.git
```
Adding to path in Matlab can be done in more than one way: 
- Per session:
https://ch.mathworks.com/help/matlab/ref/addpath.html
- Permanently at startup:
https://ch.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html

Next, get the Bioinformatics Toolbox for Matlab, more details regarding add-on management can be found here:
https://ch.mathworks.com/help/matlab/matlab_env/get-add-ons.html

Finally download AlloDy either manually from this repo or by cloning it

For the smoothest experience, AlloDy also uses visual molecular dynamics (VMD) to transform .xtc files to .dcd and to predict secondary structure of a protein. You can download it for free from here: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD 

Don't forget to add VMD to your path variable too!

You can run AlloDy without having VMD installed if you provide .dcd trajectories instead of .xtc and if you provide *helices* variable in the input scripts (*input_md2path.m* and *input_kldiv.m*).

## Instructions:

All the input parameters needed for the analysis are provided in *input_md2path.m* for **md2path** and *input_kldiv.m* for **kldiv**.
Copy the corresponding input script to your local directory and set up the input parameters. All parameters are explained in the input scripts.

Create a PDB file with separate chains for receptor, ligand, and G-protein/effector site, and define the chains in the input script.

NOTE: all input residues should be in the same pdb numbering as the input PDB

The script takes in either .xtc or .dcd trajectory files, and every trajectory file needs to be in its own directory (runN, where N is the number of the run). Prepare your directory in the following manner:

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

Run *md2pathMain.m* or *kldivMain.m*. All the output can be found in a folder called **md2path** created in your defined directory in your input script.
To run pathway calculation for individual trajectory clusters, setup and run *pca_1b_ligandPathCalcClusters.m* after the PCA and dihedrals sections in *md2pathmain.m*

Typical runtime for calculations: 
- md2path: ~ 1 hour from start to finish on a typical PC and for around 15000-20000 frames of trajectory of a 300-400 residue GPCR
- kldiv: a few minutes assuming you only calculate 1st order KL divergences

## Acknowledgements:
- Mutual information statistical filtering: McClendon et al., J Chem Theory Comput (2009)
- Using Kullback−Leibler divergences to compare conformational ensembles: McClendon et al., J Chem Theory Comput (2012)
- Clustering mutual information into allosteric pathways: Bhattacharya and Vaidehi, Biophy J (2014); Nivedha et al., Mol Pharmacol (2018)
- MDToolbox for Matlab: Matsunaga, and Sugita, J Chem Phys (2018)


## Contact: mahdi.hijazi@epfl.ch (2022)
##

