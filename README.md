# Collapse Lineage Tree

CLTree is a tool to annotate the phylogenetic tree by lineage and
measure their differences in topology by Shannon entropy.

## Program and Tasks

All tasks are run by the program `cltree` with tasks. The available tasks are:

- run: the main program, annotate the tree (newick format) by the lineage of genomes,
  and statistics. It can be emitted.
- cache: A tool to convert the NCBI taxonomy database dump file (download from:  
  [NCBI Taxonomy dump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
  to an input lineage string file
- search: Query the lineage of genomes from the NCBI database and Lineage files.
- leaf: Obtain the leaf name list of a phylogenetic tree (in Newick form).

## Installation

### Compile with CMake

#### Preparation

- cmake >= 3.0
- compiler supporting C++11 standard
- require library: libz, nlohmann-json

#### Compiling

1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make install (_option_)

### Run Programs in Docker

Docker allows users preforming programs on both Windows and Linux/MacOS.
You can download docker free and reference [docker document](https://docs.docker.com/install/)
to install it. After installing docker, basic usages for CVTree are:

1. Build/download docker image: `docker build -t="cltree-img" .`
   or `docker pull ghzuo/cltree`. In this step, an image with cltree
   programs will obtain. Here option "-t" sets the image name. After building
   the image, you can delete the dangling images for build by `docker image prune`.
2. Start container from the image in:
   `docker run --rm -it -v $PWD/example:/root/data cltree-img`
   In this step, you will enter the cltree container, and the "example" folder
   in the host will be mounted on the "data" folder in the container. Change the path
   to the data folder, and run `cltree`. You will get the result for genomes
   in the "list" file. You can change the path "\$PWD/example" to your data directory.
3. Exit and stop container: `exit` in docker terminal.
4. Run cvtree in docker by one command in example folder:
   `docker run --rm -v $PWD:/root/data cltree-img cltree`
5. More usage for docker can reference [docker document](https://docs.docker.com/).

## Run Programs with Example

If this is the first time you use CLTree package, please go to the
"example" folder. Download the taxdump.tar.gz file from
[NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).
Run the collapse command to get an annotated phylogenetic
tree and monophyly status by:

    ../build/bin/cltree

More detail of the command usages can be obtained by `-h` option. To speedup
the process, you can get the database image the data by `cltree cache` command
at first, and with option (`-N`) to output the newick format with NHX metadata
for displaying on [iTOL](https://itol.embl.de/).

## TODO

1. package database by sqlite
2. access lineage by accession number
3. UI interface program by Electron

## Reference

- Guanghong Zuo (2022) CLTree: Annotate Phylogenetic Tree by Lineage and
  Measure their Consistency based on Shannon Entropy. In preparation.
- Guanghong Zuo, Bailin Hao (2015) CVTree3 web server for
  whole-genome-based and alignment-free prokaryotic phylogeny and
  taxonomy, Genomics Proteomics & Bioinformatics, 13: 321-331

## License

This software is free for non-commercial use. For commercial use,
a software agreement is required.
