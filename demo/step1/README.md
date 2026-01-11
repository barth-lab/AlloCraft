# Step1 demo README:
This is the first step of the process, where we analyze molecular dynamics simulations to extract a wealth of information about our simulations. The relevant outputs for step2 are the allosteric hotspots and structural perturbations between states.

## Input:
Sample input file, input pdb, and trajectories are given to test run md2pathMain.m and kldivMain.m
A sample residue table for the pdb is provided in the database folder. Other reference structures will be downloaded when running md2pathMain.m or kldivMain.m

## Output: 
Expected output for given short trajectories is in md2pathdev folder in `1-D2_DA_WT`. Note that the provided trajectories are too short to output anything meaningful!!!
This also means that pathway calculation will fail starting from the sample trajectories, since entropies would not converge

## Running simulations:
gmx_input_files.zip contains all the needed files to run full length simulations on GPU machines
Simulations can also be run on CPU with proper changes to gpu_gromacs.sh

## Running the analysis:
The analysis can be executed from the interactive MATLAB GUI.

1.  **Open MATLAB.**
2.  **Navigate to the demo directory.** In the "Current Folder" toolbar, paste the absolute path to the `1-D2_DA_WT` directory and press Enter.
    ```matlab
    /path/to/AlloDy_Analysis_Repo/demo/step1/1-D2_DA_WT
    ```
3.  **Run the scripts.** In the "Command Window" (at the `>>` prompt), type the following commands. The semicolon at the end of each command suppresses verbose output and allows them to run sequentially.
    ```matlab
    input_md2path;
    md2pathMain;
    ```

## Folder description:
`1-D2_DA_WT`: Contains DD2R structure bound to dopamine (DA), which was obtained starting by PDB code 6VMS followed by docking of DA

`2-D2_BRC_WT`: Contains DD2R structure bound to bromocriptine (BRC), which was based on PDB code 6VMS

`database`: The database folder for AlloCraft. Database contains reference PDBs downloaded by AlloCraft, residue tables for generic GPCR numbering, and FASTA files containing sequences of studied proteins

`gmx_input_files.zip`: Input files for running MD simulations of DA and BRC bound DD2R using gromacs

