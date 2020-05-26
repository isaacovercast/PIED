![panda pied python](./img/pied_small.png)

# PIED
A birth/death tree with abundances. Abundance can evolve either as a BM process
with random fission at speciation events, or the rate of change of abundance
can evolve as BM.

## Default CLI args
The default CLI will parse a handful of universally useful arguments:
* `-n`  This is the flag to create a new params file
* `-p`  The flag to specify a params file for a new run of the app
* `-s`  How many simulations to run (if it's a simulation based model)
* `-c`  How many cores to spin up with the ipyparallel backend
* `-f`  Force the operation (overwrite anything that already exists)
* `-v`  Print out more progress info
* `-q`  Don't print anythin to standard out
* `-d`  Turn on debug mode to log debug info to a file
* `-V`  Print version info and exit

Long form arguments:

* `--ipcluster <cluster_id>`    Pass in a cluster ID for a running up ipcluster
    parser.add_argument("--ipcluster", metavar="ipcluster", dest="ipcluster",

## Usage
Create a params file:

    PIED -n wat

Edit the params file to update params. Run 10 simulations:

    PIED -p params-wat.txt -s 10
