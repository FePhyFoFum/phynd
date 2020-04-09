# phynd
phynd the conflict in your alignments


## General usage

The scripts here are relatively simple and easy to run. But in order to take full advantage of them, you will want to ensure that these are installed:
 - [iqtree](http://www.iqtree.org) - for running the tree analyses
 - R and ggplot within R - for plotting
 - [bp](https://github.com/FePhyFoFum/gophy) - for conducting conflict analyses

## Installation

The package is just one script, however, it relies on the three packages above. Iqtree just needs to be in your PATH (so that when you type `iqtree` in the terminal it runs). The same with R if you want to plot (also ensure the `ggplot` library is installed). `bp` is the only odd one. For this, you need to have go installed.

Use your package manager (`sudo apt install golang-go` or what is relevant) or download from the [website](https://golang.org). Then run `go get github.com/FePhyFoFum/gophy` and then `go build github.com/FePhyFoFum/gophy/bp/bp.go`. This will generate the `bp` program that you can put in your PATH (probably `sudo cp bp /usr/local/bin`).

To get phynd, you can clone the repository with `git clone https://github.com/FePhyFoFum/phynd.git` or just download the archive. 

## Running

You can run a basic analysis with the command

`./src/phynd -s SEQFILE -w 1000 -i 1000 -t 2`
and
`./src/phynd -s SEQFILE -w 1000 -i 1000 -p -t 2`
if you want to plot. 

## Examples

![Example 1](examples/t1files/t1.fa.interval_plotdata.png "Example 1")
