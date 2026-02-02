# GenKIC Dataset Creation

In order to use these scripts you need to install `PyRosetta`.

This can be done following one of these examples:
* [pyrosetta-installer](https://pypi.org/project/pyrosetta-installer/)
* [pyrosetta direct install](https://www.pyrosetta.org/downloads)

------------------------------------------------------------------------------
`insolution_genkic.py` arguments:
```
p.add_argument('-s', '--size', nargs="+", type=int, default=[6, 7, 8, 9], help="List of number of residue macrocycles you'd like to generate. (e.g. --size 6 7)")
p.add_argument('-n', '--nstruct', type=int, default=10_000, help="Number of structures to generate per size.")
p.add_argument('--debug', action="store_true", help="Dump PDBs for testing, unmute Rosetta, print helpful trace messages")
p.add_argument('--nofilter', action="store_true", help="Dont apply strict filter on hbonds")
p.add_argument('--sample-root', action="store_true", help="If generating samples in-solution then this should be set as your root residue is generally not an anchor, which means it can sample freely.")
p.add_argument("--time-test", action="store_true", help="This is only for a time comparison between XML and PyRosetta as XML writes to disk for every structure, but here we do not. This can be done for time comparisons. Dont set as it will make the process slower.")
p.add_argument("--empty-pose-test", action="store_true", help="Run an empty pose test, since this is more representative of the XML script")
p.add_argument("--alanine", action="store_true", help="Make our pose all Alanines, this is important if you want your conformations to be all L-alpha-amino acids, if not passed then pose will be all Glycines. Which allow heterochiral conformations."
```
------------------------------------------------------------------------------
Example Run:
```
python insolution_genkic.py --size 6 7 --nstruct 1000 --sample-root --alanine
```
* `--size 6 7`: will create a 6mer and 7mer ensembles
* `--nstruct`: Generate ensembles with 1000 conformations
* `--sample-root`: Apply randomization to the root residue, as GenKIC will not sample the root.
* `--alanine`: Our pose will be poly-alanine, if not set then they are poly-glycine
------------------------------------------------------------------------------

If you are looking for other cyclic peptide & design based scripts / functions in PyRosetta please checkout my public repo: ![UsefulRosettaScripts](https://github.com/PowersPope/UsefulRosettaScripts)
