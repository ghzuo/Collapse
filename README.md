# Collapse Lineage Tree

CLTree is a tool to annotate the phylogenetic tree by lineage and
measure their differences in topology by Shannon entropy.

## Program and Tasks

All tasks are run by the program `cltree` with tasks. The available tasks are:

- run: the main program, annotate the tree (newick format) by the lineage of genomes,
  and statistics. It can be emitted.
- cache: A tool to convert the NCBI taxonomy database dump file (download from:  
  [NCBI Taxonomy dump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz))
  to an input lineage string file
- search: Query the lineage of genomes from the NCBI database and Lineage files.
- leaf: Obtain the leaf name list of a phylogenetic tree (in Newick form).

## Installation

### Compile with CMake

#### Preparation

- cmake >= 3.10
- compiler supporting C++11 standard
- require library: libz, nlohmann-json

#### Compiling

1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make install (_option_)

### Run Programs in Singularity

Singularity is a container technology developed by Lawrence Berkeley National Laboratory specifically for high-performance computing (HPC) scenarios. It employs virtualization entirely based on portability, offering a more lightweight solution with faster deployment. Currently, Singularity is widely adopted by various high-performance computing centers worldwide.

1. Install singularity/apptainer, e.g., `sudo apt-get install singularity`
2. Use the cltree.sif in script directly, Or Compile your CLTree singularity container by script/cltree.def with command:
   `singularity  build --fakeroot cltree.sif cltree.def`
3. More usage for singularity can reference [singularity document](https://sylabs.io/docs/).

### Run Programs in Docker

We recommend utilizing Singularity for containerized execution of cltree. However, if you are accustomed to Docker, a Dockerfile is also provided in the script directory for compatibility,
or pull from Docker Hub by `docker pull ghzuo/cltree`.

## Run Programs with Example

If this is the first time you use CLTree package, please go to the
"example" folder. Download the taxdump.tar.gz file from
[NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).
Run the collapse command to get an annotated phylogenetic
tree and monophyly status by:

    ../build/bin/cltree

More detail of the command usages can be obtained by `-h` option. To speedup
the process, you can get the database image the data by `cltree cache` command
at first, and with option (`-I`) to output the newick tree and annotate files
for handling tree on [iTOL](https://itol.embl.de/).

## Reference

- Guanghong Zuo (2025) CLTree: Annotating, Rooting and Evaluating Phylogenetic
  Tree based on Lineage of Genomes. In preparation.
- Guanghong Zuo, Bailin Hao (2015) CVTree3 web server for
  whole-genome-based and alignment-free prokaryotic phylogeny and
  taxonomy, Genomics Proteomics & Bioinformatics, 13: 321-331

## License

This software is free for non-commercial use. For commercial use,
a software agreement is required.
