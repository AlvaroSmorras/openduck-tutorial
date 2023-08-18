# Welcome to the Openduck tutorial

Openduck is an opensource version of Dynamic Undocking, a particular implementation of steered molecular dynamics (SMD) ment to assess the robustness of protein-ligand complexes through the work needed to bring the main hydrogen bond interaction to the quasi-bound state.

## Getting started

### Installing openduck
To get started, we will install openduck on an anaconda environment (you can do it in whichever python environment you prefer).

```{bash}
$ git clone git@github.com:CBDD/openduck.git
$ cd openduck
$ conda env create -f environment.yaml
$ conda activate openduck
$ python setup.py install
$ cd ..
```

### Downloading the tutorial data

Once we have openduck up and running, we can start the tutorial by getting the necessary files to run it (including this file).

```{bash}
$ git clone git@github.com:AlvaroSmorras/openduck-tutorial.git
$ cd openduck-tutorial
```
To run Dynamic Undocking we only need a protein receptor and one or more ligands. In these tutorials, we will use ligands and proteins available in the pdb and configuration files in yaml format.

## System preparation

In this tutorial we will work with the Cyclin Dependant Kinase 2 (CDK2) and an inhibitor (SCQ), both from the [2R3K](https://www.ebi.ac.uk/pdbe/entry/pdb/2r3k) structure from the PDB. As seen in the image below, the kinase inhibitor is in a region called the 'hinge', between the beta-sheets and alpha-helixes domains. The inhibitor is interacting with the hinge loop through two hydrogen bonds (HB). Both HB are being performed with the leucine 83 backbone. Through the tutorial we will learn how to prepare and launch Dynamic Undocking to obtain the $W_{QB}$ of the inhibitor SCQ bound to CDK2.  

<p float='middle'>
<img src="./imgs/2r3k_pdb.png" width="50%" height="50%">
<img src="./imgs/ligand_HB.png" width="45%" height="45%">
</p>


### 1 Chunking

Dynamic Undocking's main descriptor, the quasi-bond work ($W_{QB}$), is a local measurement of the stability of the protein-ligand complex. This lets us reduce the receptor to the minimum residues to preserve the ligand's environment and to quicken the simulations. This proces we call chunking.

Chunking is integrated in all of the openduck preparation ( openMM-full-protocol, openMM-prepare and amber-prepare ) however, it also has a standalone subcommand. While chunking can be done within the pipeline, we recommend doing separatedly beforehand in order to check the chunk representation of the receptor and to ensure it is the most appropiate.

An *adecuated chunk* has to fulfill the following conditions:

    - All residues directly interacting with the ligand
    - No artificial gaps in the pocket close to the hydrogen bond (water shielding is one of the most
      important ligands have on WQB
    - Available ligand exit pathway (only relevant for enclosed pockets)

We will use one of the two HB between the ligand and receptor mentioned above to create the chunk around. The tutorial has been set to employ configuration files (.yaml), however, a more *classical* flag-based execution is also possible. You can find template input yamls for every subprocess in the openduck files
To create the receptor chunk we only need the receptor and ligand files, the desired interaction we will use to launch DUck and a cutoff radius.

In our case, everything is specified in the configuration file *chunk_input.yaml*

``` {yaml}
# chunk_input.yaml
# Main arguments
interaction : A_LEU_83_N
receptor_pdb : 2r3k_receptor.pdb
ligand_mol : 2r3k_lig.mol
output : 2r3k_chunk.pdb

# Chunking arguments
cutoff : 11
ignore_buffers : False
```

To launch it, we need to activate the openduck conda environment, enter the directory and execute the chunking protocol we specified in the yaml.
```
$ conda activate openduck
$ cd 1_Chunking
$ openduck chunk -y chunk_input.yaml
```

You can open the chunked receptor with your prefered visualization program to check if it has the conditions we mentioned above. As the receptor is being 'cut' to reduce the atoms, each segment needs to be capped. Check that all the segments are properly capped in the resulting receptor.
If the receptor does not fulfil the conditions of an *adecuated chunk*, you can adjust the cutoff threshold at will.

<p align='center'>
<img src="./imgs/2r3k_chunk_10A.png" width="60%" height="60%">
</p>

## 2 Parametrization

Now that we have a chunk we can proceed to parametrize the ligand, chunk and solvation. In the openduck executable this step is separated depending on the executable one wants to use afterwards, either Amber or openMM. However, the parametrization is done equally in both executions. Both *openmm-prepare* and *amber-prepare* have incorporated the chunking step, where you can use the appropiate parameters found during the previous step. Alternatively you can use the already chunked receptor but, remember checking the new interaction definition as during the chunking, the receptor residues might change numbering and chain.

For the preparation we have a plethora of options regarding forcefields, the periodic box parameters and other execution options such as hydrogen mass repartitioning (HMR). To know more about the options you can run the *openduck amber-prepare* or the *openduck openmm-prepare* with the help flag.

The preparation has the following configuration file:
```
# amber-prep_input_single_mol.yaml
# Main arguments
interaction : _LEU_31_N
receptor_pdb : 2r3k_chunk.pdb
ligand_mol : 2r3k_lig.mol

# Chunk
do_chunk : False

# Preparation
small_molecule_forcefield : gaff2
protein_forcefield : amber14-all
water_model : TIP3P
ionic_strength : 0.1
solvent_buffer_distance : 10
HMR : True

# Production arguments for amber queue and inputs
smd_cycles : 10
wqb_threshold : 6
queue_template : Slurm
```

We have chosen to prepare the solvation box with TIP3P waters, in a square box with 10A of buffer distance between the limits of the protein and the end of the box and a ionic strength of 0.1M. The chunk and ligand will be parametrized using the amber14SB and GAFF2 respectively. 
As we are going to use Amber for the DUck execution we will need to specify additional parameters. *smd_cycles* defines the amount of iterations will the protocol run for, iterations resulting in a $W_{QB}$ lower than *Wqb_threshold* will halt the execution. *HMR* will be performed on the topology to allow a 4 fs timestep. Finally, a *slurm* queue file with the DUck commands for Amber is generated. You can define new templates in the openduck files for your HPC facilities, if none (or *local*) is specified, the DUck protocol in Amber will be writen into a bash script for local execution.

```
$ cd 2_Parametrization
$ openduck amber-prepare -y amber-prep_input_single_mol.yaml
```

We now have the directory filled with different files, from the Amber input files (\*.in, dist_duck.rst & dist_md.rst), the topology & initial coordinates (HMR_system_complex.prmtop & system_complex.inpcrd) to the queue file (*duck_queue.q*). This will be all the necessary files to launch the production,but first lets have a look at the solvated system we will simulate. 

<p align='center'>
<img src="./imgs/2r3k_chunk_solvated.png" width="60%" height="60%">
</p>

### 2a Parameterizing multiple ligands

DUck was initially designed as a post-docking filter for high-througput virtual screening (HTVS) campaigns. As such, the single protein-ligand parametrization explained in the [parametrization](#2-parametrization) section is not practical when the amount of ligands to assess grow. For this reason, there is a *batch* execution for the preparation, which takes a multiligand sdf instead of a single mol and paralelizes the parametrization through the specified amount of *threads*.

The files generated will be the same as for the single ligand, albeit separated in subfolders named LIG_target_1, LIG_target_2 [..] LIG_target_n.

```
$ cd 2a_Multiple_Ligands
$ openduck amber-prepare -y amber-prep_input_multiple-ligands.yaml
```

## 3 Production

Once the system is prepared, we only need to run the simulations in you prefered machine.
The Amber commands are setup to use pmemd.cuda, which uses GPU, but in openMM we have the option to employ either CPU or GPU. Take into account that the CPU execution will be much slower. For the purpose of this tutorial, the DUck results have been precomputed and are stored in the *3_Production* directory.

The DUck pipeline has the following structure:

<p align='center'>
<img src="./imgs/duck_pipeline.png">
</p>

 
The ligand is free to explore different conformations during the equilibration and MD phase, while the receptor is restrained. The two SMD steps at different temperatures bring the specified hydrogen bond from 2.5$\AA$ to 5.0$\AA$ at a constant speed of $5\AA/\mu s$. Each SMD simulation is stored in a directory named DUCK_n or DUCK_325K_n depending on the temperature and production cycle is it. The force and work values extracted from it are stored in the duck.dat file of such directories. Production iterations are stopped prematurely if the $W_{QB}$ does not reach the specified *wqb_threshold* after each SMD.

## 4 Analysis

There are various ways of analyzing the results. The most quick and straightforward is using the min $W_{QB}$, as it represents the most probable pathway being the least resistant. The $W_{QB}$ of each replica is usually reflected at the last step of the simulation. However, this is not always the case so it needs to be recovered from the maximum work after the minima at ~3$\AA$. A more thorough approach is obtaining the quasi-bond free energy $\Delta G_{QB}$ by applying the Jarzynski equality (JE).

<p align='center'>
<img src="./imgs/wqb_dqb_schema.png" width="80%">
</p>


The JE is a particular case of the fluctuation-dissipation theorem. In a non-equilibrium simulation, transition between a two states (in our case between the bound and quasi-bound state) dissipates energy, turning it into heat (i.e. friction).This friction increases the faster the process goes, and is 0 when the process is infinitely slow. In this particular case, the work (*W*) needed for the process equates the free energy.

$$\Delta F = W - W_{diss} $$ 

Through the JE we relate the quasi-bond free energy with the Boltzmann average of the works obtained during out SMD simulations.Normally, the more production cycles, the more accurate is the $\Delta G_{QB}$. This however, can be solved by artificially subsampling the simulations through bootstrapping, for example.

$$e^{{-\Delta F}/{kT}} = \overline{e^{-W/kT}}$$

$$\Delta G_{QB} = -kTln( M^{-1} \sum_{i=1}^{M}{exp({-W_{i}}/{kT})}$$

### 4a Analysis of a single ligand

With the simulations we produced in the [previous section](#3-production) we will analyze the structural stability of our protein-ligand complex. First, with a simple command, we will calculate and plot the $W_{QB}$. We will specify the flag *-d* to change the output to yield an average and standart deviation (SD) of the $W_{QB}$ and the *--plot* to obtain a graphical representation of work during the simulations. 

```{bash}
$ cd 4_Analysis
$ openduck report -d avg --plot
$ eog wqb_plot.png 
```
|System	|WQB	|Average	|SD|
|-|-|-|-|
|.|	7.268099999999999	|8.186807272727272	|0.409612311140689|
As you can see, the $W_{QB}$ is 7.268 and the simulations present a very small hysteresis. This is good and shows the ligand has a stable binding mode (low *SD*) and a robust interaction with the receptor (high $W_{QB}$). 

<p align='center'>
<img src="./4_Analysis/wqb_plot.png" width="70%">
</p>

 Take notice, that the average in this table does not correspond with the $\Delta G_{QB}$, as it is a arithmetic average and not an Boltzmann average. To calculate the jarzynski free energy we can invoke the report command with *jarzynski* in the *-d* flag.

 ```{bash}
 $ openduck report -d jarzynski
 ```
|System	|Jarzynski	|Jarzynski_SD	|Jarzynski_SEM|
|-|-|-|-|
|.|	8.0316067705202	|0.10848616866492665|0.017371689901702|
<p align='center'>
<img src="./4_Analysis/bootstraped_WQB_plot.png" width="70%">
</p>

### 4b Highthroughput analysis of multiple ligands


