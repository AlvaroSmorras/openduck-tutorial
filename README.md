# Welcome to the Openduck tutorial

Openduck is an opensource version of Dynamic Undocking, a particular implementation of steered molecular dynamics ment to assess the robustness of protein-ligand complexes through the work needed to bring the main hydrogen bond interaction to the quasi-bound state.

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

In this tutorial we will work with the Cyclin Dependant Kinase 2 (CDK2), more specifically the [2R3K](https://www.ebi.ac.uk/pdbe/entry/pdb/2r3k) structure from the PDB.

![cdk2_receptor](./imgs/2r3k_pdb.png)

### 1_Chunking

Dynamic Undocking's main descriptor, the quasi-bond work ($W_{QB}$), is a local measurement of the stability of the protein-ligand complex. This lets us reduce the receptor to the minimum residues to preserve the ligand's environment and to quicken the simulations. This proces we call chunking.

Chunking is integrated in all of the openduck preparation ( openMM-full-protocol, openMM-prepare and amber-prepare ) however, it also has a standalone subcommand. While chunking can be done within the pipeline, we recommend doing separatedly beforehand in order to check the chunk representation of the receptor and to ensure it is the most appropiate.

The optimal chunk has to fulfill the following conditions:

    - All residues directly interacting with the ligand
    - No artificial gaps in the pocket close to the hydrogen bond (water shielding is one of the most important ligands have on $W_{QB}$)
    - Available ligand exit pathway (only relevant for enclosed pockets)


```
$ conda activate openduck
$ cd 1_Chunking
$ openduck chunk -y chunk_input.yaml
```

![cdk2_chunk](./imgs/2r3k_chunk.png)