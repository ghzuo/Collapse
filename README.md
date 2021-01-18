# Collapse Lineage Tree

A tool to annotate phylogenetic tree by lineage and measure their differences in toplogy by Shannon entropy.

## Program and Tasks

All tasks are run by the program `cltree` with tasks. The available tasks are:

- run: the main program, annotate the tree (newick format) by lineage,
  and statistics.
- index: A tool to convert the NCBI taxonomy database dump file (download from:  
  [NCBI Taxonomy dump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
  to an input lineage string file
- query: Query the lineage of genome names/taxids from the database.
- leaf: Obtain the leafs name list of a phylogenetic tree (in Newick form).

## Installation

### Compile with CMake

#### Preparation

- cmake >= 3.0
- compiler supporting C++11 standard
- require library: libz

#### Compiling

1. unzip the package file and change into it
2. mkdir build and change into it
3. cmake .. or add some options you wanted
4. make
5. make install (_option_)

### Run Programms in Docker

Docker allows users run programs on both Windows and Linux/MacOS.
You can download docker free and reference [docker document](https://docs.docker.com/install/)
to install it. After install docker, basic usages for CVTree are:

1. Build/download docker image: `docker build -t="cltree" .`
   or `docker pull ghzuo/cltree`. In this step, a image with cvtree
   programs will obtained. Here option "-t" set the image name. After build
   image, you can delete the dangling images for build by `docker image prune`.
2. Start container from image:
   `docker run --rm -it -v $PWD/example:/root/data cltree`
   In this step, you will enter the cvtree container, and the "example" folder
   of this project will be find in the "data" folder. Change path to the data folder,
   and run `cltree`. You will get the result for eight genomes in the "list"
   file. You can change the path "\$PWD/example" to your own data directory.
3. Exit and stop container: `exit` in docker terminal.
4. Run cvtree in docker by one command:
   `docker run --rm -v $PWD:/data -w /data cltree cltree`
5. More usage for docker can reference [docker document](https://docs.docker.com/).

## Run Programs with Example

If this is the first time you use Collapse package, please go to the
"example" folder. Download the taxdump.tar.gz file from
[NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) to
this dirtory. Run the collapse command to get an annotated phylogenetic
tree and monophyly status by:

    ../build/collapse

More detail of the command usage can be obtaion by `-h` option. To speedup
the process, you can get the database image the data by `getdb` command at
first.

## TODO

1. detail usage
2. package database by sqlite
3. access lineage by accession number
4. UI interface program by Electron

## Reference

- Guanghong Zuo, Bailin Hao (2015) CVTree3 web server for
  whole-genome-based and alignment-free prokaryotic phylogeny and
  taxonomy, Genomics Proteomics & Bioinformatics, 13: 321-331

## License

This software is free for non-commercial use. For commercial use,
a software agreement is required.
