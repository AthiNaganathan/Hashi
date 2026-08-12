# _Hashi_: Bridging Statistical Model Derived 1D Microstate Encodings and Protein 3D Structural Ensembles

## About

_Hashi_ is a pipeline to rapidly and efficiently generate realistic structural ensembles from the outputs of the structure-based Wako-Saitô-Muñoz Eaton (WSME) statistical mechanical model of protein folding. This approach relies on integrating the block WSME model outputs – strings of zeros and ones describing the conformational status of every residue over thousands or millions of microstates each assigned a statistical weight derived from physically grounded energy-entropy terms, and free energy profiles - with the RANCH module of the EOM (ensemble optimization method) from the ATSAS software suite, providing three-dimensional views of the structural ensembles within the model framework.

In this README, we cover the installation and usage requirements for _Hashi_, with specific examples.

## Installation Instructions

### Software instructions

Following version of softwares are required -
- MATLAB v24.2.0.2773142 (R2024b) Update 2
    - Bioinformatics Toolbox
    - Optimization Toolbox
    - Parallel Computing Toolbox
    - Signal Processing Toolbox
    - Statistics and Machine Learning Toolbox
- ATSAS v4.1.0
-  PyMOL(TM) 3.1.6.1 (Only used for visualization of final outputs)

> **WARNING -** There are some known inconsistencies of this pipeline with more recent version of ATSAS. The user is recommended to install specifically the version mentioned above.

#### Note for Linux/MacOS users -

Explicit initializing of environment variables is required in Linux/MacOS prior to calling terminal from MATLAB. 

If running _Hashi_ on MacOS, in `scripts\init.m`, delete line 12 and update environment variables ATSAS and PATH in lines 13 and 14 (commented guide available in the code, an error is thrown up if this is not set up properly).

If running _Hashi_ on Linux, in `scripts\init.m`, delete line 16 and update environment variables ATSAS, ATSAS_LICENSE, PATH and LD_LIBRARY_PATH in lines 17-20 (commented guide available in the code, an error is thrown up if this is not set up properly).

#### Note for generating STRIDE output -

The bWSME model requires a STRIDE output for calculating block assignments. As mentioned in the subsequent section, this output can be generated from https://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py. However, should the server not be available, Linux/MacOS users (or Windows users with WSL) can generate the same by downloading and running STRIDE locally from https://github.com/MDAnalysis/stride. While it is advisable to refer to the official documentation for installation specifications, a brief summary of the same would look as -

```
git clone https://github.com/MDAnalysis/stride.git
cd stride/src
make
./stride your_protein.pdb
```

And save the output to a text file as `struct.txt`.

## Usage

### Overview

The overall workflow is split into two steps - (1) Running bWSME calculations to generate 1D free energy profiles and ensembles, and (2) passing this information to RANCH to generate 3D structural ensembles.

> On cloning the repo, with `\scripts\` as the working directory, run `init.m` to generate the empty directories `\data\` and `\results\` which are required for successful execution of the code.

> All scripts are to be run with `\scripts\` being the working directory. All the scripts that are intended to be directly executed by the user are present in this directory exclusively.

A general run requires three input files. For a protein named `pname`, these are -
1. `pname.pdb` - PDB structure of the protein
2. `pname.cif` - CIF structure of the protein. Must have the same chain ID and atoms as the PDB
3. `struct.txt` - The formatted STRIDE output (refer note in previous section) should be saved as a text file which will be read by the code to generate blocks

> **NOTE -** The user must ensure good quality input .pdb and .cif files for the program, with no missing atoms and consistent residue indices. Incomplete or corrupted input files have been known to loop the RANCH algorithm indefinitely, or give rise to other errors.
>
> Further, presence of HETATOM records in the input files also causes errors in the algorithm, so the user is advices to remove these atoms beforehand.

With these files in `\data\`, `main_WSME.m` is called. The outputs of `main_WSME.m` is stored in `\data\WSME_outputs\`. Two specific output files, `pname_pepval.mat` and `pname_BlockDet.mat` are automatically copied into `\data\`, as these two files (along with `pname.pdb` and `pname.cif`) are the four input files required for `main_RANCH.m`. The user is recommended to properly read and execute the two examples given in the next section, to understand the flow of the program, before running it on their data.

> **NOTE :** If you already have `pname_pepval.mat` and `pname_BlockDet.mat` in `\data\` (say, from a previous run), you can run `main_RANCH.m` directly.

Following this, the input parameters can be modified in `main_RANCH.m` (explained in examples), and the script can be executed. The outputs can be found in `\results\`. 
- `\results\microstates\` contains all ranked microstate information listed sequentially (as an intermediate step in executing RANCH).
- `\results\pools\` is the intended final output folder, with conformers sorted into user-defined pools.

> **A Note on terminologies used -**
> 1. A "Pool" is defined as a range of macrostates (say, between 19 to 25 structured blocks). Pools are user-defined.
> 2. A "Macrostate" is defined on the basis of the number of structured blocks. All conformations with the same number of structured blocks belong to the same macrostate.
> 3. A "Microstate" is defined by the location of structured/unstructured blocks along the protein sequence. The output of the bWSME model is the thermodynamic ranking of all possible microstates of a protein within the model framework.
> 4. A "Confromation" is the unit 3D structure. A microstate can have different conformations that differ in the physical location of the atoms in the structured/unstructured regions. Conformations are generated for each microstate by the RANCH algorithm, and are physically reasonable with no steric clashes. Note that only the Calpha atoms are generated for the unstructured regions by the RANCH algorithm.

Alongside the structures, `\results\pools\` also contains the information of the microstates that constitute that pool. Each row corresponds to one of the microstates in the pool, and the format of the columns is as follows -

>| Column number | Value interpretation |
>| :--- | :--- |
>| **1** | Number of structured blocks in the microstate |
>| **2** | The overall probability of this microstate against all microstates in the scope of the bWSME model |
>| **3** | The starting block index of the first structured region |
>| **4** | The ending block index of the first structured region |
>| **5** | The starting block index for the second structured region (set to 0 for SSA) |
>| **6** | The ending block index for the second structured region (set to 0 for SSA) |
>| **7**  | Set to 1 for SSA, 2 for DSA, and 3 for DSAw/L
>| **8** | The rank of this microstate when ranked amongst all microstates considered in that run (across all pools). The microstate with rank 1 has the highest probability of occuring amongst all microstate across all pools defined in that specific run.

This can be used as a reference for further analysis, if necessary.

> **WARNING :** Every run of `main_WSME.m` clears the `\data\WSME_outputs\` folder, and every run of `main_RANCH.m` clears the `\results\` folder. Make sure to save these outputs, should you want to. The scripts will prompt you prior to the deletion of these files, at the start of every run.

## Example 1: Villin

The following example uses data from `\examples\Villin\` to demonstrate the typical workflow. Villin is a short protein (35 resideus), and the following example considers only one pool. **Example 2: Pertactin** considers a larger protein with multiple pools.

### 1.1 Setting Up -

If not already done, with `\scripts\` as the working directory, run `\scripts\init.m` to generate the empty directories `\data\` and `\results\`. 

Copy the three files (`Villin.pdb`, `Villin.cif` and `struct.txt`) from `\examples\Villin\` into `\data\`. The directory should now have these three files along with the empty directory `\data\WSME_outputs\`.

### 1.2 Running main_WSME.m -

Open `\scripts\main_WSME.m` and edit the following parameters -
1. Set `pname` to the protein name (line 13), `BlockSize` to the size of blocks to be considered in the bWSME model (line 14), and `is_long_protein` to zero (line 15). 
2. Set the input parameters for `cmapCalcElecBlock` and `FesCalc_Block` (their documentation can be found in the function source codes in `\scripts\pipelines`) in the parameters between lines 16 and 23. 

For the ease of the user, these values default to the current Villin example, given as -
```matlab
% Villin
pname  = 'Villin';
BlockSize = 1;
is_long_protein = 0;
struct_file = 'struct.txt';
pH = 7;
srcutoff = 5.0;
ene = -87.5/1000;
DS = -14.5/1000;
DCp = -0.36/1000;
T = 310;
IS = 0.15;
```
The parameter `is_long_protein` determines whether `FesCalc_Block()` (for short proteins, such as Villin) or `FesCalc_Block_gen()` (for longer proteins, such as Pertactin in Example 2) is called for determining the free energy profile. If `is_long_protein` is set to 1, definition of additional parameters `disr` and `ppos` is required, as shown in the next example.

>  **NOTE :** It is empirically observed that for proteins longer than ~200 residues, setting `is_long_protein` to 1 leads to significant performance boost. Exact performance depends on the value of other input parameters as well.

3. Then, execute `main_WSME.m`. When prompted, type `y` at the command line and hit `Enter`.

You will now see output files of the same in `\data\WSME_outputs\`, alongside two new files `Villin_pepval.mat` and `Villin_BlockDet.mat` that are automatically copied into `\data\`. These files are required for running `main_RANCH.m`.

Calling function `disp_blocks(pname)` from the MATLAB command line displays the 1D Macrostate landscape, which can serve as a guide for assigning pools of interest. For the current Villin example, the output from calling function `disp_blocks('Villin')` looks as follows -

![Villin - Visualization of 1D landscape](/images/Villin_fig1.png)


### 1.3 Running main_RANCH.m -

Open `\scripts\main_RANCH.m` and edit the following parameters -
1. Set `pname` to the protein name (line 16).
2. Set `macrostate_pools` to the pools you wish to visualize (line 17).
3. Set `no_of_microstates` to the number of (statistically most significant) microstates to be considered per pool (line 18).
4. Set `no_of_conformations` to the number of conformations to be generated per microstate (line 19).

For the ease of the user, these values default to the current Villin example, given as -
```matlab
% Example 1 : Villin
pname = 'Villin'; % Protein name
macrostate_pools = [25 34]; % Split 1D landscape into pools wherein to visualize ensemble
no_of_microstates = 20; % No. of microstates considered per pool of macrostates
no_of_conformations = 5; % No. of conformations to be generated per microstate
```

Once pools have been assigned, calling function `disp_blocks(pname, macrostate_pools)` displayes the assigned pools, for clarification. For this example, the output is as follows -

![Villin - Assignment of Pools](/images/Villin_fig2.png)

5. Once satisfied with the pools assignment, execute `main_RANCH.m`. When prompted, type `y` at the command line and hit `Enter`.

The output pool and microstate information can be found in `\results\pools\`. The structures can then be visualized as required using softwares such as PyMol, VMD, etc.

The naming of the output conformers follows the pattern `<no of structured blocks>_<SSA(1) or DSA(2) or DSAw/L(3)>_<rank across all pools>_<order of generation by RANCH>.cif`

### 1.4 Outputs -

> **NOTE -** Different runs will generate different conformations due to random seed initialization by RANCH. 

The outputs in `\results\pools\pool_1\` for the above example are visualized below. These are possible conformations of the top 20 most probable microstates with 25-34 structured residues.

![Villin - Outputs in pool 1](/images/Villin_fig3.png)

The user is of course encouraged to utilize their own visualization pipelines. The above is just for a quick visual.


## Example 2: Pertactin

### 2.1 Setting Up -

If not already done, with `\scripts\` as the working directory, run `\scripts\init.m` to generate the empty directories `\data\` and `\results\`.

Copy the three files (`Pertactin.pdb`, `Pertactin.cif`, and `struct.txt`) from `\examples\Pertactin\` into `\data\`. The directory should now have these three files along with the empty directory `\data\WSME_outputs\`.

### 2.2 Running main_WSME.m -

Open `\scripts\main_WSME.m` and edit the following parameters -
1. Set protein name and edit input parameters as described in the previous example. For a long proteins (empirically, >150 residues) such as Pertactin, it is recommended to set `is_long_protein` to 1.

For the ease of the user, the values for this run are available commented from the lines **(25-36)**. You can comment the values of the previous example (lines **12-23**), and uncomment the lines for this example. Once uncommented, they appear as follows -
```matlab
% Pertactin
pname = 'Pertactin';
BlockSize = 5;
is_long_protein = 1;
struct_file = 'struct.txt';
pH = 7;
srcutoff = 5.0;
ene = -79/1000;
DS = -14.5/1000;
DCp = -0.36/1000;
T = 298;
IS = 0.05; 
```

3. Then, execute `main_WSME.m`. When prompted, type `y` at the command line and hit `Enter`.

You will now see output files of the same in `\data\WSME_outputs\`, alongside two new files `Pertactin_pepval.mat` and `Pertactin_BlockDet.mat` in `\data\`. The visaulization of the 1D landscape can be done by calling `disp_blocks('Pertactin')` as in the previous example.

### 2.3 Running main_RANCH.m -

Open `\scripts\main_RANCH.m` and edit the input parameters. The meanings of the various input parameters remain identical to the previous example.

For the ease of the user, the values for this run are available commented from the lines **(21-25)**. You can comment the values of the previous example (lines **15-19**), and uncomment the lines for this example. Once uncommented, they appear as follows -
```matlab
% Example 2 : Pertactin
pname = 'Pertactin'; % Protein name
macrostate_pools = [68 76; 107 111; 128 137]; % Split 1D landscape into pools wherein to visualize ensemble
no_of_microstates = 5; % No. of microstates considered per pool of macrostates
no_of_conformations = 5; % No. of conformations to be generated per microstate
```

> **WARNING -** It is a known issue for this protein that when RANCH is asked to model very short double Glycine-containing sequences between two larger structured blocks (in the DSAw/L approximation), the RANCH program gets stuck for arbitrarily long times. This behaviour so far has only occured in this protein (for certain microstates only in pools like 11-15 or 86-95, not covered in this example) out of the 30+ proteins tested. Resolving this would require more detailed examination of the implementation of the RANCH algorithm. 
>
>It is possible, though unlikely, that similar behaviour could be observed on other proteins where specific sequences in a specific ordered/disordered schemes cause RANCH to get stuck indefinitely. Currently, the variable `timeoutSec` in `line 13` (set by default to 5 minutes, since RANCH runs typically take <2 minutes) controls how long RANCH is allowed to run for on any particular microstate before forcibly being stopped. The microstates where the conformations did not converge are stored in `results/skipped_mis.mat` (in standard format, except the 8th column of each entry points to the pool from where that microstate was taken), and a summary of the skipped microstates are displayed in the MATLAB console after the execution of the run. The user may choose to increase `timeoutSec` if they feel their device might be slower, but in our testing, if RANCH did not converge in ~2 min, it did not converge even after multiple hours.The `plot_fails` function is useful in visualizing the microsates that were skipped due to timeout.

Once pools have been assigned, calling function `disp_blocks(pname, macrostate_pools)` displayes the assigned pools, for verification. For this example, the output is as follows -

![Pertactin - Assignment of Pools](/images/Pertactin_fig2.png)

5. Once satisfied with the pools assignment, execute `main_RANCH.m`. When prompted, type `y` at the command line and hit `Enter`.

The output pool and microstate information can be found in `\results\pools\`. The structures can then be visualized as required using softwares such as PyMol, VMD, etc.

The naming of the output conformers follows the pattern `<no of structured blocks>_<SSA(1) or DSA(2) or DSAw/L(3)>_<rank across all pools>_<order of generation by RANCH>.cif`

### 2.4 Outputs -

> **NOTE -** Different runs will generate different conformations due to random seed initialization by RANCH. 

The outputs in `\results\pools\pool_3\` for the above example are visualized below. These are possible conformations of the top 5 most probable microstates with 128-137 structured residues, corresponding to the folded ensemble.

![Villin - Outputs in pool 1](/images/Pertactin_fig5.png)

The outputs in `\results\pools\pool_2\` for the above example are visualized below. These are possible conformations of the top 5 most probable microstates with 107-111 structured residues, corresponding to the intermediate-like ensemble.

![Villin - Outputs in pool 1](/images/Pertactin_fig4.png)

The outputs in `\results\pools\pool_1\` for the above example are visualized below. These are possible conformations of the top 5 most probable microstates with 68-76 structured residues, corresponding to the intermediate ensemble.

![Villin - Outputs in pool 1](/images/Pertactin_fig3.png)

The script for aligning+coloring structures used by the author is provided in `\scripts\` as `align_structures.py`. The residues to be aligned can be modified in lines 19 and 20 (default value is set to align the C terminus residues 366-536 for the current Pertactin example). When called from PyMol, it loads, aligns and colours all the .cif files present in the working directory. To use it, copy it into any of the `\results\pools\pool_n\` directories, and, with the same being the working directory in PyMol,  call it from PyMol using `run align_structures.py` from the PyMol command line. The output is saved as `aligned_session.pse` in the same directory.

The user is of course encouraged to utilize their own visualization pipelines. The above is just for a quick visual.

## References

If using _Hashi_, please cite the following two articles:

1. Soundhararajan Gopi, Akashnathan Aranganathan & Athi N. Naganathan (2019). Thermodynamics and Folding Landscapes of Large Proteins from a Statistical Mechanical Model. Curr. Res. Struct. Biol., 1, 6-12.

2. Athi N. Naganathan, Rahul Dani, Soundhararajan Gopi, Akashnathan Aranganathan & Abhishek Narayan (2021). Folding Intermediates, Heterogeneous Native Ensembles and Protein Function. J. Mol. Biol., 433, 167325.

Additional reading and references:

1. Sathvik Anantakrishnan & Athi N. Naganathan (2023). Thermodynamic Architecture and Conformational Plasticity of GPCRs. Nat. Commun., 14, 128.

2. Athi N. Naganathan & Adithi Kannan (2021). A Hierarchy of Coupling Free Energies Underlie the Thermodynamic and Functional Architecture of Protein Structures. Curr. Res. Struct. Biol., 3, 257-267.

3. Soundhararajan Gopi, Animesh Singh, Swaathiratna Suresh, Suvadip Paul, Sayan Ranu & Athi N. Naganathan (2017). Toward a Quantitative Description of Microscopic Pathway Heterogeneity in Protein Folding. Phys. Chem. Chem. Phys., 19, 20891 - 20903.

4. Nandakumar Rajasekaran, Swaathiratna Suresh, Soundhararajan Gopi, Karthik Raman & Athi N. Naganathan (2017). A General Mechanism for the Propagation of Mutational Effects in Proteins. Biochemistry, 56, 294–305.

5. Athi N. Naganathan (2012). Predictions from an Ising-Like Statistical Mechanical Model on the Dynamic and Thermodynamic Effects of Protein Surface Electrostatics. J. Chem. Theory Comput., 8, 4646-4656.
