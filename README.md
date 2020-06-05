![panda pied python](./img/pied_small.png)

# PIED
A birth/death process with abundances and genetic diversity at the tips. 
Abundance can evolve either as a BM process with random fission at speciation
events, or the rate of change (r) of abundance can evolve as BM in which case
abundance (n) changes through time (dt) via `n * exp(r\*dt)`. Speciation rate
can also shift at branching events in the manner of ClaDS. 

## Installation

* Install [conda](https://docs.conda.io/en/latest/miniconda.html)  for python3
* `conda create -n PIED python=3.7`
* `conda activate PIED`
* `conda install -c conda-forge -c iovercast pied`

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
    0.1                  ## [12] [lambda_sigma]: Rate at which speciation rate changes if ClaDS is True.
    0.1                  ## [13] [alpha]: Rate shift if ClaDS is True
    500                  ## [14] [sequence_length]: Length of the genomic region simulated, in base pairs.
    1e-05                ## [15] [mutation_rate]: Mutation rate per base per generation
    10                   ## [16] [sample_size]: Number of samples to draw for calculating genetic diversity

Run 10 simulations:

    PIED -p params-wat.txt -s 10

Run 10 simulations on 10 cores in parallel:

    PIED -p params-wat.txt -s 10 -c 10

## Output
Results are written to `default_PIED/wat-SIMOUT.csv`. Not generally human
readable the results file contains the parameters used to generate each
simulation, as well as the calculated extinction rate (as a fraction of
birth events), data for each tip, including abundance, genetic diversity,
growth rate, and speciation rate, and finally the dated tree in newick form.
Field names are as follows (along with 1 example simulation):

    birth_rate stop_criterion ntaxa time process ClaDS abundance_mean abundance_sigma growth_rate_mean growth_rate_sigma lambda_sigma alpha sequence_length mutation_rate sample_size ext_rate data tree
    1 time 20 4.0 abundance False 50000 0.1 0.0 0.01 0.1 0.1 500 1e-05 10 0.0 r13:108:0.003644444444444445:-0.0015646353607632964:0.9820914460240241,r12:123:0.005422222222222223:-0.0020408591149106933:0.9837358548596271,r11:1590:0.08671111111111121:0.0010199726182759323:0.9803538911026921,r10:57:0.0007111111111111111:0.0014860494178264385:1.014698708129729,r9:553:0.02231111111111111:0.0010251014167980722:1.0296906741616945,r8:253:0.006311111111111113:-0.001134422411248719:0.9831400703641143,r7:1590:0.04182222222222222:-6.416658964779895e-05:1.0471523763965862,r6:8808:0.1276:0.003977195684297366:0.973644150250009,r5:590:0.025555555555555547:-0.0017646281295530285:1.0028299834839174,r4:55:0.0013333333333333335:-0.0025752837037666693:1.0141098013462215,r3:4446:0.17804444444444473:-0.002410783979211202:0.9593471649972787,r2:21782:0.5767999999999978:0.0003952310179508536:0.9804683207023276,r1:10413:0.34182222222222136:0.008606087644055268:1.0406933431412086,r0:401:0.01648888888888889:-0.009038279539788011:0.8548065098437879 ((((r13:0.116531,r12:0.116531)0:0.206095,r11:0.322626)0:0.823284,((r10:0.217932,r9:0.217932)0:0.0945876,r8:0.31252)0:0.83339)0:2.86134,(((r7:0.568332,r6:0.568332)0:0.394463,(r5:0.829775,r4:0.829775)0:0.13302)0:1.296,((r3:0.36856,r2:0.36856)0:1.80003,(r1:1.83062,r0:1.83062)0:0.337969)0:0.0902007)0:1.74846);

The data for each simulation can be parsed in python like this:

    simfile = "/path/to/default_PIED/wat-SIMOUT.csv"
    sim_df = pd.read_csv(simfile, header=0, sep=" ")
    sims = []
    # There's probably a "fancier" way to do this, but this gets the job done
    for rec in df["data"]:
        # split the records for each species, separated by ','
        dat = rec.split(",")
        # Create a dictionary with species id as the key and the value is another
        # dictionary mapping 'abundance', 'pi', 'r', & 'lambda_' to their respective
        # values per species.
        dat = {x:{"abundance":int(y), "pi":float(z), "r":float(aa), "lambda_":float(bb)} for x, y, z, aa, bb in map(lambda x: x.split(":"), dat)}
        # append this dictionary to the sims list
        sims.append(dat)

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
