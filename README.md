# PIED
A birth/death tree with abundances. Abundance can evolve either as a BM process
with random fission at speciation events, or the rate of change of abundance
can evolve as BM.

## Where the work is done
Inside PIED.__main__ around line 159 is a call to `data.simulate()`. This
will call the workhorse function of the primary class, and so PIED.Core
should have a `simulate()` function which does the work and accepts and
handles well the following arguments:
* nsims     The number of simulations to run 
* ipyclient The ipyparallel client to use
* quiet     Don't print anything to standard out
* verbose   Print more stuff to standard out
* force     Overwrite previous simulations if necessary

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
* Clone this repo
* Replace PIED with your chosen new package name:

    egrep -lRZ 'PIED' . | xargs -0 -l sed -i -e 's/PIED/foo/g'

* If you choose, you can rename the `Core` class to something more meaningful

    egrep -lRZ 'Core' . | xargs -0 -l sed -i -e 's/Core/Bar/g'

* You are probably 95% there, but there will probably be some things to clean
up. Better than doing it all from scratch.
