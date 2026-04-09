# Hashi : Bridging Binary Microstate Encodings and Protein 3D Structural Ensembles

## About

Description of Hashi

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

### OS-specific dependencies

## Usage

### Overview

The overall workflow is split into two steps - (1) Running bWSME calculations to generate 1D free energy profiles and ensembles, and (2) passing this information to RANCH to generate 3D structural ensembles.
The overall workflow is split into two steps - (1) Running bWSME calculations to generate 1D free energy profiles and ensembles, and (2) passing this information to RANCH to generate 3D structural ensembles.

> On cloning the repo, with `\scripts\` as the working directory, run `\scripts\init.m` to generate the empty directories `\data\` and `\results\` which are required for successful execution of the code.

> All scripts are to be run with `\scripts\` being the working directory. All the scripts that are intended to be directly executed by the user are present in this directory exclusively.

A general run requires three input files. For a protein named `pname`, these are -
1. `pname.pdb` - PDB structure of the protein
2. `pname.cif` - CIF structure of the protein. Must have the same chain ID and atoms as the PDB
3. `struct.txt` - The formatted STRIDE output - the output from feeding the structure into the STRIDE websever (https://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py) should be saved as a txt file which will be read by the code to generate blocks

With these files in `\data\`, `main_WSME.m` is called. The outputs of `main_WSME.m` is stored in `\data\WSME_outputs\`. Two specific output files, `pname_pepval.mat` and `pname_BlockDet.mat` are automatically copied into `\data\`, as these two files (along with `pname.pdb` and `pname.cif`) are the four input files required for `main_RANCH.m`. The user is recommended to properly read and execute the two examples given in the next section, to understand the flow of the program, before running it on their data.

> **NOTE :** If you already have `pname_pepval.mat` and `pname_BlockDet.mat` (say, from a previous run), you can run `main_RANCH.m` directly, as shown in **Example 2 : Pertactin**.

Following this, the input parameters can be modified in `main_RANCH.m` (explained in examples), and the script can be executed. The outputs can be found in `\results\`. 
- `\results\microstates\` contains all microstate information listed sequentially (as an intermediate step in executing RANCH).
- `\results\pools\` is the intended final output folder, with conformers sorted into user-defined pools.

> **A Note on terminologies used -**
> 1. A "Pool" is defined as a range of macrostates, and are user-defined. All conformations having between 19-25 structured blocks/residues (for example) can be placed in one pool.
> 2. A "Macrostate" is defined on the basis of the number of structured blocks/residues. All conformations with the same number of structured blocks/residues belong to the same macrostate.
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

### Setting Up -

If not already done, with `\scripts\` as the working directory, run `\scripts\init.m` to generate the empty directories `\data\` and `\results\`. 

Copy the three files (`Villin.pdb`, `Villin.cif` and `struct.txt`) from `\examples\Villin\` into `\data\`. The directory should now have these three files along with the empty directory `\data\WSME_outputs\`.

### Running main_WSME.m -

Open `\scripts\main_WSME.m` and edit the following parameters -
1. Set `pname` to the protein name (line 9). Here, set
```matlab
% Set protein name
pname = 'Villin';
```
2. Set the input parameters for `cmapCalcElecBlock` and `FesCalc_Block` (their documentation can be found in the function source codes in `\scripts\pipelines`) (lines 33 and 35). Here, set
```matlab
% Inputs = protein name, stride output file, pH, srcutoff, BlockSize
cmapCalcElecBlock(pname, 'struct.txt', 7, 5.0, 1)
% Inputs = protein name, stride output file, ene, DS, DCp, T, IS
FesCalc_Block(pname, 'struct.txt', -98/1000, -14.5/1000, -0.3579/1000, 310, 0.1)
```

3. Then, execute `main_WSME.m`. When prompted, type `y` at the command line and hit `Enter`.

You will now see output files of the same in `\data\WSME_outputs\`, alongside two new files `Villin_pepval.mat` and `Villin_BlockDet.mat` in `\data\`.

Calling function `disp_blocks(pname)` from the MATLAB command line displays the 1D Macrostate landscape, which can serve as a guide for assigning pools of interest. For the current Villin example, the output from calling function `disp_blocks('Villin')` looks as follows -

![Villin - Visualization of 1D landscape](/Villin_fig1.png)


### Running main_RANCH.m -

Open `\scripts\main_RANCH.m` and edit the following parameters -
1. Set `pname` to the protein name (line 10).
2. Set `macrostate_pools` to the pools you wish to visualize (line 11).
3. Set `no_of_microstates` to the number of (statistically most significant) microstates to be considered per pool (line 12).
4. Set `no_of_conformations` to the number of conformations to be generated per microstate (line 13).

For the ease of the user, these values default to the current Villin example, given as -
```matlab
% Example 1 : Villin
pname = 'Villin'; % Protein name
macrostate_pools = [25 34]; % Split 1D landscape into pools wherein to visualize ensemble
no_of_microstates = 20; % No. of microstates considered per pool of macrostates
no_of_conformations = 5; % No. of conformations to be generated per microstate
```

Once pools have been assigned, calling function `disp_blocks(pname, macrostate_pools)` displayes the assigned pools, for verification. For this example, the output is as follows -

![Villin - Assignment of Pools](/Villin_fig2.png)

5. Once satisfied with the pools assignment, execute `main_RANCH.m`. When prompted, type `y` at the command line and hit `Enter`.

The output pool and microstate information can be found in `\results\pools\`. The structures can then be visualized as required using softwares such as PyMol, VMD, etc.

The naming of the output conformers follows the pattern `<no of structured blocks>_<SSA(1) or DSA(2) or DSAw/L(3)>_<rank across all pools>_<order of generation by RANCH>.cif`

## Example 2: Pertactin

### Setting Up -

If not already done, with `\scripts\` as the working directory, run `\scripts\init.m` to generate the empty directories `\data\` and `\results\`.

The directory should now have these three files along with the empty directory `\data\WSME_outputs\`.

>**NOTE :** In this example, we omit the generation of the `Pertactin_pepval.mat` and `Pertactin_BlockDet.dat` files, choosing instead to provide these outputs directly in the `\examples\Pertactin` directory. This is because Pertactin is 538 residues long, and it requires ~400GB of runtime memory to compute the same. The `struct.txt` file is provided, however, should the user have access to the computing power required to generate the same.

Copy the four files (`Pertactin.pdb`, `Pertactin.cif`, `Pertactin_pepval.mat` and `Pertactin_BlockDet.dat`) from `\examples\Pertactin\` into `\data\`.

### Running main_RANCH.m -

Open `\scripts\main_RANCH.m` and edit the input parameters. The meanings of the various input parameters remain identical to the previous example.

For the ease of the user, the values for this run are available commented from the lines **(15-19)**. You can comment the values of the previous example (lines **9-13**), and uncomment the lines for this example. Once uncommented, they appear as follows -
```matlab
% Example 2 : Pertactin
pname = 'Pertactin'; % Protein name
macrostate_pools = [68 76; 107 111; 128 137]; % Split 1D landscape into pools wherein to visualize ensemble
no_of_microstates = 5; % No. of microstates considered per pool of macrostates
no_of_conformations = 5; % No. of conformations to be generated per microstate
```

> **WARNING -** It is a known issue that increasing `no_of_conformations` to values beyond ~5 (sometimes) causes the program to halt indefinitely, especially for large proteins as in this example. This is because RANCH does not timeout, it searches till it finds a solution, and for unlucky initial conditions, the algorithm gets stuck in a loop. The fix for this is currently not available, and the user is recommended to be conservative while assigning the input parameters, and to retry execution, should the program get stuck.

Once pools have been assigned, calling function `disp_blocks(pname, macrostate_pools)` displayes the assigned pools, for verification. For this example, the output is as follows -

![Pertactin - Assignment of Pools](/Pertactin_fig2.png)

5. Once satisfied with the pools assignment, execute `main_RANCH.m`. When prompted, type `y` at the command line and hit `Enter`.

The output pool and microstate information can be found in `\results\pools\`. The structures can then be visualized as required using softwares such as PyMol, VMD, etc.

The naming of the output conformers follows the pattern `<no of structured blocks>_<SSA(1) or DSA(2) or DSAw/L(3)>_<rank across all pools>_<order of generation by RANCH>.cif`

## Notes

## References
