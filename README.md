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
simulation, as well as the observed numbers of tips, the observed simulation
time, and the calculated extinction rate (as a fraction of birth events),
data for each tip, including abundance, genetic diversity, growth rate,
and speciation rate, and finally the dated tree in newick form. Field names
are as follows (along with 1 example simulation):

    birth_rate stop_criterion ntaxa time process ClaDS abundance_mean abundance_sigma growth_rate_mean growth_rate_sigma lambda_sigma alpha sequence_length mutation_rate sample_size obs_ntaxa obs_time ext_rate data tree
    1 time 20 4.0 abundance False 50000 0.1 0.0 0.01 0.1 0.1 500 1e-05 10 23 4.095569867467191 0.0 r22:6866:0.3388888888888878:-0.0017392248027676867:0.9929405172387488,r21:3409:0.06493333333333326:0.0022979523708463396:1.044889056623618,r20:1361:0.04755555555555551:-0.0008304076216114711:1.028788037227819,r19:2153:0.047066666666666625:-0.0010936878297818674:0.9892171283717283,r18:2071:0.13764444444444443:-0.03602844187649056:1.2236414272712746,r17:1794:0.033111111111111105:0.010652702902916795:0.8841139779859748,r16:1600:0.04968888888888883:0.00022248012819702791:1.007759644044947,r15:206:0.016355555555555557:-0.0025001157661682927:1.0119988921069114,r14:121:0.0032:0.0028114408498106655:1.0150306188449219,r13:7176:0.1911555555555557:-0.0002989927924176867:1.0080428685670018,r12:57:0.0026222222222222224:0.0022534864488625425:1.0310102680219881,r11:233:0.015377777777777773:-0.0016532087475619054:0.9874506470688181,r10:2650:0.040133333333333306:0.0012761459270739734:0.9590447837930026,r9:1072:0.0538666666666666:0.0017936880382974293:1.047772584251838,r8:41:0.0010666666666666667:-0.00155993687616825:1.0010437225970643,r7:571:0.01613333333333333:-0.0011302498187538252:1.027807449254878,r6:131:0.005466666666666667:-0.0034603340148914993:0.9977304525473011,r5:2790:0.04586666666666664:-0.003898338884046415:0.9882789905775415,r4:972:0.022755555555555553:0.0009140839738056091:0.9665357418700201,r3:9040:0.27791111111111194:0.00047583092183889905:0.9715740561602696,r2:668:0.04315555555555552:0.0005252530117281093:0.9894320705866041,r1:4840:0.12311111111111138:-0.0002025660774967209:1.0244034979593213,r0:531:0.013777777777777774:0.001325242277879697:0.9853079841162011 ((r22:1.06953,(r21:0.410235,(r20:0.348699,r19:0.348699)0:0.0615366)0:0.659294)0:3.02604,((r18:3.28678,r17:3.28678)0:0.150255,((((r16:0.356944,((r15:0.1582,r14:0.1582)0:0.0917129,r13:0.249913)0:0.107031)0:0.388738,(r12:0.562355,(r11:0.169675,r10:0.169675)0:0.39268)0:0.183328)0:0.894806,(((r9:0.351402,(r8:0.170165,r7:0.170165)0:0.181237)0:0.395484,r6:0.746886)0:0.39907,(r5:0.7366,r4:0.7366)0:0.409356)0:0.494533)0:0.216775,(r3:0.860269,((r2:0.447829,r1:0.447829)0:0.328436,r0:0.776265)0:0.0840038)0:0.996995)0:1.57977)0:0.658533);

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

* `--ipcluster <cluster_id>`    Pass in the cluster ID of a running ipcluster instance
