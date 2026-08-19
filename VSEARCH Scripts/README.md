This folder contains the two Linux scripts used to clean Nanopore sequencing outputs including; deprelication, filtering, chimera removal.

The primary script, eDNA_pipeline.sh, contains the full pipeline from start to final VSEARCH. However, some steps at the beginning are commented out as these had been performed already - remove comments as necessary

The second script, temp_vsearch.sh, just focuses on running a VSEARCH command - this script was created to run VSEARCH without re-running the entire pipeline.

Reference FASTAs cited within scripts may differ - as such double check before use and edit file names to use the desired FASTA/file.
