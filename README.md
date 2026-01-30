# CyclicMPNN - A cyclic peptide specialized fine-tuned version of ProteinMPNN
<p align="center">
  <img src="https://drive.google.com/uc?export=view&id=11DKu21MMIl0KgKobTAUcPOg7BHXchi8H" width="800">
</p>

Read [CyclicMPNN paper]().

This repository is a derivative of [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) (Dauparas et al., 2022)

We generated in silico data and fine-tuned ProteinMPNN, resulting in an overall improvement in cyclic peptide sequence
design performance.

## What is different from ProteinMPNN?
- Fine-tuned on cyclic peptide dataset (*in silico* & PDB; see paper)
- Fine-tuning training procedure described in Methods
- No architectural changes
- Same inference pipeline

## Abstract
Cyclic peptides are a promising class of therapeutics due to their attractive drug qualities such
as increased structural stability, cell permeability, and resistance to proteolytic degradation. With
recent advancements in cyclic peptide backbone generation models like CyclicCAE and RFPeptide,
generating cyclic peptide backbones can be done more rapidly compared to traditional algorithm
or physics-based approaches. However, designing energetically favorable cyclic peptide sequences
to fit generated backbones using only canonical amino acids is nontrivial. We fine-tuned the state-
of-the-art deep learning model for protein sequence design, ProteinMPNN, using a combination
of x-ray crystal structures from the Protein Data Bank and *in silico* generated cyclic peptides. Our
approach surpasses ProteinMPNN in cyclic peptide sequence design, producing energetically stable
sequences with a higher success rate of folding into the generated cyclic peptide backbones. We show
that CyclicMPNN can be used as a motif-inpainting strategy and in *de novo* sequence design tasks.
We propose that CyclicMPNN will enable the rapid design of energetically stable cyclic peptide
sequences, increasing the success rate of therapeutic cyclic peptide development.

-----------------------------------------------------------------------------------------------------
To run CyclicMPNN clone this github repo and install Python>=3.0, PyTorch, Numpy. 

Full cyclic peptide model weights: `cyclicmpnn_weights/cyclicmpnn_48_010.pt`

Helper scripts to generate poly-alanine *in silico* ensembles: `helper_scripts/insolution_genkic.py`

Code organization:
* `protein_mpnn_run.py` - the main script to initialialize and run the model.
* `protein_mpnn_utils.py` - utility functions for the main script.
* `examples/` - simple code examples.
* `inputs/` - input PDB files for examples
* `outputs/` - outputs from examples
* `colab_notebooks/` - Google Colab examples
* `training/` - code and data to retrain the model
* `cyclicmpnn_weights/` - specific CyclicMPNN weights
* `helper_scripts/` - CyclicMPNN *in silico* generation PyRosetta scripts and utils

## CyclicMPNN Citation
```
@article{,
  title={CyclicMPNN: Stable Cyclic Peptide Sequence Generation},
  author={Powers, Andrew Carl and Janthana, Yanapat and Hosseinzade, Parisa},
  journal={},
  volume={},
  number={},  
  pages={},
  year={2026},
  publisher={}
}
```
-----------------------------------------------------------------------------------------------------
For example to make a conda environment to run CyclicMPNN:
* `conda create --name cyclicmpnn` - this creates conda environment called `cyclicmpnn`
* `source activate cyclicmpnn` - this activate environment
* `conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch` - install pytorch following steps from https://pytorch.org/
-----------------------------------------------------------------------------------------------------
Output example:
```
>size10_ALA_bb_008430, score=1.3696
LGEGGPGTLR
>size10_ALA_bb_003908, score=1.5825
WVTKSDNSEY
>size10_ALA_bb_005886, score=1.4018
RYAHPTSAES
>size10_ALA_bb_006418, score=1.4044
RTGPADSDNA
```
* `score` - average over residues that were designed negative log probability of sampled amino acids

## Original ProteinMPNN Documentation

This repository uses the unmodified ProteinMPNN codebase.

For full documentation of scripts, flags, and examples, see:
[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)

## Original ProteinMPNN Citation
```
@article{dauparas2022robust,
  title={Robust deep learning--based protein sequence design using ProteinMPNN},
  author={Dauparas, Justas and Anishchenko, Ivan and Bennett, Nathaniel and Bai, Hua and Ragotte, Robert J and Milles, Lukas F and Wicky, Basile IM and Courbet, Alexis and de Haas, Rob J and Bethel, Neville and others},
  journal={Science},
  volume={378},
  number={6615},  
  pages={49--56},
  year={2022},
  publisher={American Association for the Advancement of Science}
}
```
-----------------------------------------------------------------------------------------------------
