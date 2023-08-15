# kb_jorg release notes
=========================================


0.5.0
-----
* Relaxed MIRA naming stringency
* Major Jorg re-write in python for improved KBase logging
* Added multi-threading parameter to MIRA

0.4.1
-----
* Reorganized main jorg function
* Decluttered variable usage throughout
* Start Jorg python re-write

0.4.0
-----
* Added minimap2 and bwa as selectable read mappers
* Option to keep all Iterations folder fna files

0.3.0
-----
* Added auto-renaming for long sequence IDs
* Internalized Jorg binaries

0.2.0
-----
* Base image overhaul
* Allow users to input multiple assemblies

0.1.1
-----
* Lower kmer minimum size

0.1.0
-----
* Release to beta

0.0.1
-----
* Under development

0.0.0
-----
* Module created by kb-sdk init

ToDo
_____
* reinstate runtime cap
* fix tests so can be run in sequence, required intermittent file deletion
* reinstate mira and mirabait logs
* review Jorg repo for more workflow additions
* save circularized contig information in assembly object
* single-end read mode
* make Jorg tutorial
* write a manuscript
