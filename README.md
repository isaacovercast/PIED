![panda pied python](./img/pied_small.png)

# PIED
A birth/death process with abundances and genetic diversity at the tips. 
Abundance can evolve either as a BM process with random fission at speciation
events, or the rate of change (r) of abundance can evolve as BM in which case
abundance (n) changes through time (dt) via `n * exp(r)\*\*dt`. Speciation rate
can also shift at branching events in the manner of ClaDS. 

## Usage
Create a params file:

    PIED -n wat

Look at the params and edit them if you wish:

    ------- PIED params file (v.0.0.2)----------------------------------------------
    wat                  ## [0] [simulation_name]: The name of this simulation scenario
    ./default_PIED       ## [1] [project_dir]: Where to save files
    1                    ## [2] [birth_rate]: Speciation rate
    taxa                 ## [3] [stop_criterion]: Whether to stop on ntaxa or time
    20                   ## [4] [ntaxa]: Number of taxa to simulate if stop is `ntaxa`
    4                    ## [5] [time]: Amount of time to simulate if stop is `time`
    abundance            ## [6] [process]: Whether to evolve `abundance` or growth `rate` via BM
    False                ## [7] [ClaDS]: Whether to allow speciation rates to change along the branches a la ClaDS
    50000                ## [8] [abundance_mean]: Ancestral abundance at time 0
    0.1                  ## [9] [abundance_sigma]: Rate at which abundance changes if process is `abundance`
    0                    ## [10] [growth_rate_mean]: Ancestral population growth rate at time 0.
    0.01                 ## [11] [growth_rate_sigma]: Rate at which growth rate changes if process is `rate`
    1                    ## [12] [lambda_mean]: Ancestral speciation rate at time 0.
    0.1                  ## [13] [lambda_sigma]: Rate at which speciation rate changes if ClaDS is True.
    0.1                  ## [14] [alpha]: Rate shift if ClaDS is True
    1e-05                ## [15] [mutation_rate]: Mutation rate per base per generation
    10                   ## [16] [sample_size]: Number of samples to draw for calculating genetic diversity

Run 10 simulations:

    PIED -p params-wat.txt -s 10

Results are written to `default_PIED/wat-SIMOUT.csv`.

Run 10 simulations on 10 cores in parallel:

    PIED -p params-wat.txt -s 10 -c 10

## Default CLI args
The default CLI will parse a handful of universally useful arguments:
* `-n`  This is the flag to create a new params file
* `-p`  The flag to specify a params file for a new run of the app
* `-s`  How many simulations to run
* `-c`  How many cores to spin up with the ipyparallel backend
* `-f`  Force the operation (overwrite anything that already exists)
* `-v`  Print out more progress info
* `-q`  Don't print anything to standard out
* `-d`  Turn on debug mode to log debug info to a file
* `-V`  Print version info and exit

Long form arguments:

* `--ipcluster <cluster_id>`    Pass in a cluster ID for a running up ipcluster
    parser.add_argument("--ipcluster", metavar="ipcluster", dest="ipcluster",
