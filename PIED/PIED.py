import datetime
import json
import logging
import msprime
import numpy as np
import os
import string
import time
import tempfile
import toytree

from collections import OrderedDict

import PIED
from PIED.util import *

LOGGER = logging.getLogger(__name__)

class Core(object):
    """
    The PIED object

    :param str name: The name of this PIED simulation. This is used for
        creating output files.
    :param bool quiet: Don't print anything ever.
    :param bool verbose: Print more progress info.
    """

    def __init__(self, name, quiet=False, verbose=False):
        if not name:
            raise PIEDError(REQUIRE_NAME)

        ## Do some checking here to make sure the name doesn't have
        ## special characters, spaces, or path delimiters. Allow _ and -.
        ## This will raise an error immediately if there are bad chars in name.
        self._check_name(name)
        self.name = name

        self._version = PIED.__version__

        ## stores default ipcluster launch info
        self._ipcluster = {
            "cluster_id" : "",
            "profile" : "default",
            "engines" : "Local",
            "quiet" : 0,
            "timeout" : 120,
            "cores" : 0, #detect_cpus(),
            "threads" : 2,
            "pids": {},
            }

        ## the default params dict
        ## If you add a parameter to this dictionary you need
        ## to also add a short description to the PARAMS dict
        ## at the end of this file
        ##
        ## Also be sure to add it to _paramschecker so the type gets set correctly
        ## FIXME: PIED parameters need to be added here
        self.paramsdict = OrderedDict([
                       ("simulation_name", name),
                       ("project_dir", "./default_PIED"),
                       ("birth_rate", 1),
                       ("stop_criterion", "taxa"),
                       ("ntaxa", 20),
                       ("time", 4),
                       ("process", "abundance"),
                       ("ClaDS", False),
                       ("abundance_mean", 50000),
                       ("abundance_sigma", 0.1),
                       ("growth_rate_mean", 0),
                       ("growth_rate_sigma", 0.01),
                       ("lambda_mean", 1),
                       ("lambda_sigma", 0.1),
                       ("alpha", 0.1),
                       ("mutation_rate", 1e-5),
                       ("sample_size", 10)
        ])

        ## Separator to use for reading/writing files
        self._sep = " "

        ## A dictionary for storing taxon specific information. This dictionary
        ## is populated when empirical data is loaded. If no per taxon info
        ## is present then we sample from the priors as specified in the params
        ## file.
        self.taxa = {}

        ## elite hackers only internal dictionary, normally you shouldn't mess with this
        self._hackersonly = dict([
        ])

        ## The empirical msfs
        self.empirical_msfs = ""


    #########################
    ## Housekeeping functions
    #########################
    def __str__(self):
        return "<PIED.Core: {}>".format(self.paramsdict["simulation_name"])


    ## Test Core name is valid and raise if it contains any special characters
    def _check_name(self, name):
        invalid_chars = string.punctuation.replace("_", "")\
                                          .replace("-", "")+ " "
        if any(char in invalid_chars for char in name):
            raise PIEDError(BAD_PIED_NAME.format(name))


    ## This doesn't do anything, but it could be useful if you need to record
    ## per simulation data.
    def _get_simulation_outdir(self, prefix=""):
        """
        Construct an output directory for a simulation run.
        Make output directory formatted like <output_dir>/<name>-<timestamp><random 2 digit #>
        This will _mostly_ avoid output directory collisions, right?

        :param string prefix: The directory within which to create the
            simulation output directory.
        """

        dirname = prefix + self.paramsdict["simulation_name"]
        outdir = os.path.join(self.paramsdict["project_dir"],\
                              dirname\
                              + "-" + str(time.time()).replace(".", "")[-7:]\
                              + str(np.random.randint(100)))
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        return outdir


    def _paramschecker(self, param, newvalue, quiet=True):
        """
        Check and set parameters. Raises exceptions when params are set to 
        values they should not be.

        :param string param: The parameter to set.
        :param newvalue: The value of the parameter.
        :param bool quiet: Whether to print info.
        """
        ## TODO: This should actually check the values and make sure they make sense
        ## FIXME: PIED parameters need to be updated here.
        try:
            ints = ["birth_rate", "ntaxa", "abundance_mean", "sample_size"]
            floats = ["time", "abundance_sigma", "growth_rate_mean", "growth_rate_sigma",\
                        "lambda_mean", "lambda_sigma", "alpha", "mutation_rate"]
            ## Cast params to correct types
            if param == "project_dir":
                ## If it already exists then just inform the user that we'll be adding
                ## more simulations to the current project directory
                if " " in newvalue:
                    raise PIEDError("`project_dir` may not contain spaces. You put:\n{}".format(newvalue))
                self.paramsdict[param] = os.path.realpath(os.path.expanduser(newvalue))

                if not os.path.exists(self.paramsdict["project_dir"]):
                    os.mkdir(self.paramsdict["project_dir"])
            
            elif param == "stop_criterion":
                if newvalue not in ["taxa", "time"]:
                    raise PIEDError("Bad parameter: `stop_criterion` must "\
                                    + "be one of `taxa` or `time`.")
            elif param == "process":
                if newvalue not in ["abundance", "rate"]:
                    raise PIEDError("Bad parameter: `process` must be "\
                                    + "`abundance` or `rate`.")

            elif param in ints + floats:
                dtype = float
                if param in ints:
                    dtype = int
                tup = tuplecheck(newvalue, dtype=dtype)
                self.paramsdict[param] = tup
            else:
                self.paramsdict[param] = newvalue
        except Exception as inst:
            ## Do something intelligent here?
            raise


    ## Getting parameters header and parameters carves off
    ## the simulation name and the project directory
    def _get_params_header(self):
        return list(self.paramsdict.keys())[2:]


    def _get_params_values(self):
        return list(self.paramsdict.values())[2:]


    def set_param(self, param, value, quiet=True):
        """
        A convenience function for setting parameters in the API mode, which
        turns out to be a little annoying if you don't provide this. With
        the set_param method you can set parameters on the Region, the
        Metacommunity, or the LocalCommunity. Simply pass the parameter
        name and the value, and this method identifies the appropriate target
        parameter.

        :param string param: The name of the parameter to set.
        :param value: The value of the parameter to set.
        :param bool quiet: Whether to print info to the console.
        """
        try:
            self = set_params(self, param, value, quiet)
        except:
            raise PIEDError("Bad param/value {}/{}".format(param, value))


    def get_params(self):
        """
        A convenience function for getting nicely formatted params in API mode.

        :return: A string of all the params ready to be printed.
        """
        tf = tempfile.NamedTemporaryFile()
        self.write_params(outfile=tf.name, force=True)
        dat = open(tf.name).read()
        return dat


    def write_params(self, outfile=None, outdir=None, force=False):
        """
        Write out the parameters of this model to a file properly formatted as
        input for the PIED CLI. A good and simple way to share/archive 
        parameter settings for simulations. This is also the function that's
        used by __main__ to generate default params.txt files for `PIED -n`.

        :param string outfile: The name of the params file to generate. If not
            specified this will default to `params-<Region.name>.txt`.
        :param string outdir: The directory to write the params file to. If not
            specified this will default to the project_dir.
        :param bool force: Whether to overwrite if a file already exists.
        """
        if outfile is None:
            outfile = "params-"+self.paramsdict["simulation_name"]+".txt"

        ## If outdir is blank then default to writing to the project dir
        if outdir is None:
            outdir = self.paramsdict["project_dir"]
        elif not os.path.exists(outdir):
            raise PIEDError(NO_OUTDIR).format(outdir)

        outfile = os.path.join(outdir, outfile)
        ## Test if params file already exists?
        ## If not forcing, test for file and bail out if it exists
        if not force:
            if os.path.isfile(outfile):
                raise PIEDError(PARAMS_EXISTS.format(outfile))

        with open(outfile, 'w') as paramsfile:
            ## Write the header. Format to 80 columns
            header = "------- PIED params file (v.{})".format(PIED.__version__)
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            ## Whip through the current paramsdict and write out the current
            ## param value, the ordered dict index number. Also,
            ## get the short description from paramsinfo. Make it look pretty,
            ## pad nicely if at all possible.
            for key, val in self.paramsdict.items():
                if isinstance(val, tuple):
                    paramvalue = "-".join(map(str, val))
                else:
                    paramvalue = str(val)

                padding = (" "*(20-len(paramvalue)))
                paramkey = list(self.paramsdict.keys()).index(key)
                paramindex = " ## [{}] ".format(paramkey)
                LOGGER.debug("{} {} {}".format(key, val, paramindex))
                name = "[{}]: ".format(key)
                description = PARAMS[key]
                paramsfile.write("\n" + paramvalue + padding + \
                                        paramindex + name + description)

            paramsfile.write("\n")


    def save(self):
        _save_json(self)


    @staticmethod
    def load(json_path, quiet=False):
        """
        Load a json serialized object and ensure it matches to the current model format.
        """
        # expand HOME in JSON path name
        json_path = json_path.replace("~", os.path.expanduser("~"))

        # raise error if JSON not found
        if not os.path.exists(json_path):
            raise PIEDError("""
                Could not find saved model file (.json) in expected location.
                Checks in: [project_dir]/[assembly_name].json
                Checked: {}
                """.format(json_path))

        # load JSON file
        with open(json_path, 'rb') as infile:
            fullj = json.loads(infile.read(), object_hook=_tup_and_byte)

        # get name and project_dir from loaded JSON
        oldname = fullj["model"].pop("name")
        olddir = fullj["model"]["paramsdict"]["project_dir"]
        oldpath = os.path.join(olddir, os.path.splitext(oldname)[0] + ".json")

        # create a fresh new Assembly
        null = PIED.Core(oldname, quiet=True)

        # print Loading message with shortened path
        if not quiet:
            oldpath = oldpath.replace(os.path.expanduser("~"), "~")
            print("  loading Core: {}".format(oldname))
            print("  from saved path: {}".format(oldpath))

        ## FIXME: Anything more sophisticated will have to be done
        ##        custom.
        ## get the taxa. Create empty taxa dict of correct length
        #taxa = fullj["model"].pop("taxa")
        #null.taxa = taxa

        ## Set params
        oldparams = fullj["model"].pop("paramsdict")
        for param, value in oldparams.items():
            null.set_param(param, value)

        oldhackersonly = fullj["model"].pop("_hackersonly")
        null._hackersonly = oldhackersonly

        null._sep = fullj["model"].pop("_sep")

        #taxon_names = list(fullj["taxa"].keys())

        return null


    def parallel_simulate(self, ipyclient, nsims=1, quiet=False, verbose=False):
        parallel_jobs = {}
        _ipcluster = {}
        ## store ipyclient engine pids to the Core so we can
        ## hard-interrupt them later if Core is interrupted.
        ## Only stores pids of engines that aren't busy at this moment,
        ## otherwise it would block here while waiting to find their pids.
        _ipcluster["pids"] = {}
        for eid in ipyclient.ids:
            engine = ipyclient[eid]
            if not engine.outstanding:
                pid = engine.apply(os.getpid).get()
                _ipcluster["pids"][eid] = pid

        lbview = ipyclient.load_balanced_view()
        for i in range(nsims):
            ## Call do_serial sims args are: nsims, quiet, verbose
            parallel_jobs[i] = lbview.apply(serial_simulate, self, 1, True, False)

        ## Wait for all jobs to finish
        start = time.time()
        printstr = " Performing Simulations    | {} |"
        while 1:
            try:
                fin = [i.ready() for i in parallel_jobs.values()]
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                if not quiet: progressbar(len(fin), sum(fin), printstr.format(elapsed))
                time.sleep(0.1)
                if len(fin) == sum(fin):
                    break
            except KeyboardInterrupt as inst:
                print("\n    Cancelling remaining simulations.")
                break
        if not quiet: progressbar(100, 100, " Finished {} simulations in   {}\n".format(i+1, elapsed))

        faildict = {}
        param_df = pd.DataFrame()
        result_list = []
        ## Gather results
        for result in parallel_jobs:
            try:
                if not parallel_jobs[result].successful():
                    faildict[result] = parallel_jobs[result].metadata.error
                else:
                    result_list.extend(parallel_jobs[result].result())
            except Exception as inst:
                LOGGER.error("Caught a failed simulation - {}".format(inst))
                ## Don't let one bad apple spoin the bunch,
                ## so keep trying through the rest of the asyncs
        LOGGER.debug(faildict)

        return result_list
    
    
    def serial_simulate(self, nsims=1, quiet=False, verbose=False):
        import pandas as pd

        tree_list = []

        printstr = " Performing Simulations    | {} |"
        start = time.time()
        for i in range(nsims):
            try:
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                if not quiet: progressbar(nsims, i, printstr.format(elapsed))

                res = self._simulate()
                tree_list.append(res)

            except KeyboardInterrupt as inst:
                print("\n    Cancelling remaining simulations")
                break
            except Exception as inst:
                LOGGER.debug("Simulation failed: {}".format(inst))
                raise PIEDError("Failed inside serial_simulate: {}".format(inst))

        if not quiet: progressbar(100, 100, " Finished {} simulations in   {}\n".format(i+1, elapsed))

        return tree_list


    def _simulate(self, 
                alpha=0, ClaDS=False, verbose=False):
        
        # Each species has a dictionary of features, these are the default values
        #  abundance - Abundance of the species
        #  r - Rate at which abundance changes, can be negative
        #  trait - This is a random trait value that evolves by BM
        #  lambda - Per lineage speciation rate
        feature_dict = {"abundance":{"sigma":0.1, "zbar_0":50000, "log":True, "dtype":"int"},
                          "r":{"sigma":0.01, "zbar_0":0, "log":False, "dtype":"float"},
                          "trait":{"sigma":2, "zbar_0":0, "log":False, "dtype":"float"},
                          "lambda_":{"sigma":0.1, "zbar_0":self.paramsdict["birth_rate"], "log":False, "dtype":"float"}
                         }

        tre = toytree.tree()
        for fname, fdict in feature_dict.items():
            tre.treenode.add_feature(fname, fdict["zbar_0"])

        taxa_stop = self.paramsdict["ntaxa"]
        time_stop = self.paramsdict["time"]

        ext = 0
        evnts = 0
        t = 0
        while(1):

            ## Get list of extant tips
            tips = tre.treenode.get_leaves()

            # Sample time interval
            if ClaDS:
                lambs = np.array([tip.lambda_ for tip in tips])
                # Run a horse race for all lineages, smallest time sampled wins
                ts = np.random.exponential(lambs)
                idx = np.where(ts == ts.min())[0][0]
                dt = ts.min()
                sp = tips[idx]
            else:
                # This isn't strictly necessary, because if the rate shift model
                # is turned off then all tips will have equal lambda and equal
                # probability of being sampled, but keeping this here for now
                # for the hell of it.
                dt = np.random.exponential(1/(len(tips) * self.paramsdict["birth_rate"]))
                sp = np.random.choice(tips)

            t = t + dt
            evnts += 1

            c1 = sp.add_child(name=str(t)+"1", dist=0)
            c2 = sp.add_child(name=str(t)+"2", dist=0)
            try:
                # Always have at least 1 individual for each daughter species
                # Raises ValueError if sp.abundanc == 1
                abund = np.random.randint(1, sp.abundance)
            except ValueError as inst:
                # If a species of abundance 1 speciates we need to allow this
                sp.abundance = 2
                abund = 1

            for c in [c1, c2]:
                for fname, fdict in feature_dict.items():
                    c.add_feature(fname, fdict["zbar_0"])
                
            c1.add_feature("abundance", abund)
            c2.add_feature("abundance", sp.abundance-abund)

            # Update branch length, evolve features, and update abundance
            tips = tre.treenode.get_leaves()
            for x in tips:
                x.dist += dt
                if self.paramsdict["process"] == "abundance":
                    # If 'abundance' log species abundances change via Brownian motion
                    for fname, fdict in feature_dict.items():
                        x.add_feature(fname, _bm(getattr(x, fname),
                                                fdict["sigma"],
                                                dt=dt,
                                                log=fdict["log"],
                                                dtype=fdict["dtype"]))
                else:
                    for fname, fdict in feature_dict.items():
                        # Update 'lambda_', 'r' and 'trait', but skip 'abundance'
                        #import pdb; pdb.set_trace()
                        if fname == "abundance": continue
                        x.add_feature(fname, _bm(getattr(x, fname),
                                                fdict["sigma"],
                                                dt=dt,
                                                log=fdict["log"],
                                                dtype=fdict["dtype"]))
                        # Apply the population size change
                        x.abundance = int(x.abundance * (np.exp(x.r**dt)))

            # Check boundary conditions for abundance and lambda
            for x in tips[:]:
                # Check extinction and prune out extinct and any resulting hanging branches
                # If you're the only individual left, you have nobody to mate with so you die.
                if x.abundance <= 1:
                    tips.remove(x)
                    if len(tips) == 0:
                        # Can't prune to an empty tree, so simply return a new empty tree
                        tre = toytree.tree()
                    else:
                        tre.treenode.prune(tips, preserve_branch_length=True)
                    ext += 1

                # The speciation rate can't be negative
                elif x.lambda_ < 0:
                    x.lambda_ = 0

            ## This shouldn't be necessary now that we're using the treenode.prune() method
            tre = _prune(tre)

            tips = tre.treenode.get_leaves()
            # Check stopping criterion
            done = False
            if self.paramsdict["stop_criterion"]  == "taxa":
                if len(tips) >= self.paramsdict["ntaxa"]:
                    done = True
            else:
                if t >= self.paramsdict["time"]:
                    done = True
            if len(tips) == 1 and tips[0].name == '0':
                print("All lineages extinct")
                done = True
            if done:
                if verbose:
                    print("ntips {}".format(len(tips)))
                    print("time {}".format(t))
                    print("Birth events {}".format(evnts))
                    print("Extinctions (per birth) {} ({})".format(ext, ext/evnts))
                tre._coords.update()
                for i, t in enumerate(tips[::-1]):
                    t.name = "r{}".format(i)
                return tre

    
    def simulate(self, nsims=1, ipyclient=None, quiet=False, verbose=False, force=False):
        """
        Do the heavy lifting here. 

        :param int nsims: The number of PIED simulations to perform
        :param ipyparallel.Client ipyclient: If specified use this ipyparallel
            client to parallelize simulation runs. If not specified simulations
            will be run serially.
        :para bool quiet: Whether to display progress of these simulations.
        :para bool verbose: Display a bit more progress information.
        :param bool force: Whether to append to or overwrite results from
            previous simulations. Setting `force` to ``True`` will overwrite
            any previously generated simulation in the `project_dir/{name}-SIMOUT.txt`
            file.
        """
        if not quiet: print("    Generating {} simulation(s).".format(nsims))

        if not os.path.exists(self.paramsdict["project_dir"]):
            os.mkdir(self.paramsdict["project_dir"])

        if ipyclient:
            result_list = self.parallel_simulate(ipyclient, nsims=nsims, quiet=quiet, verbose=verbose)
        else:
            # Run simulations serially
            result_list = self.serial_simulate(nsims=nsims, quiet=quiet, verbose=verbose)

        self._write_df(result_list, force=force)


    ## Save the results to the output DataFrame
    def _write_df(self, result_list, force=False):

        simfile = os.path.join(self.paramsdict["project_dir"], "{}-SIMOUT.csv".format(self.name))
        ## Open output file. If force then overwrite existing, otherwise just append.
        mode = 'a'
        if force:
            ## Prevent from shooting yourself in the foot with -f
            try:
                mode = 'w'
                os.rename(simfile, simfile+".bak")
            except FileNotFoundError:
                ## If the simfile doesn't exist catch the error and move on
                pass

        with open(simfile, mode) as output:
            output.write("\n".join([x.write() for x in result_list]) + "\n")


def serial_simulate(model, nsims=1, quiet=False, verbose=False):
    import os
    LOGGER.debug("Entering sim - {} on pid {}\n{}".format(model, os.getpid(), model.paramsdict))
    res = model.serial_simulate(nsims, quiet=quiet, verbose=verbose)
    LOGGER.debug("Leaving sim - {} on pid {}".format(model, os.getpid()))
    return res

## Brownian motion function
def _bm(mean, var, dt, log=True, dtype="int"):
    ret = 0
    mean = np.float(mean)
    if dtype == "int":
        # Avoid log(1)
        if not log or mean <= 1:
            ret = np.int(np.round(np.random.normal(mean, var*dt)))
        else:
            try:
                ret = np.int(np.round(np.exp(np.random.normal(np.log(mean), var*dt))))
            except:
                import pdb; pdb.set_trace()

    elif dtype == "float":
        ret = np.random.normal(mean, var*dt)
    else:
       raise Exception("bm dtype must be 'int' or 'float'. You put: {}".format(dtype))
    return ret


def _prune(tre, verbose=False):
    ttree = tre.copy()
    tips = ttree.treenode.get_leaves()
    
    if np.any(np.array([x.height for x in tips]) > 0):
        for t in tips:
            if not np.isclose(t.height, 0):
                if verbose: print("Removing node/height {}/{}".format(t.name, t.height))
                t.delete(prevent_nondicotomic=False)
                ttree = prune(ttree)
    return ttree


##########################################
## Saving functions to dump model to json
## This is all ripped directly from ipyrad
## with minor modifications.
##########################################

class _Encoder(json.JSONEncoder):
    """
    Save JSON string with tuples embedded as described in stackoverflow
    thread. Modified here to include dictionary values as tuples.
    link: http://stackoverflow.com/questions/15721363/

    This Encoder Class is used as the 'cls' argument to json.dumps()
    """
    def encode(self, obj):
        """ function to encode json string"""
        def hint_tuples(item):
            """ embeds __tuple__ hinter in json strings """
            if isinstance(item, tuple):
                return {'__tuple__': True, 'items': item}
            if isinstance(item, list):
                return [hint_tuples(e) for e in item]
            if isinstance(item, dict):
                return {
                    key: hint_tuples(val) for key, val in item.items()
                }
            else:
                return item
        return super(_Encoder, self).encode(hint_tuples(obj))


def _default(o):
    print(o)
    # https://stackoverflow.com/questions/11942364/
    # typeerror-integer-is-not-json-serializable-when-
    # serializing-json-in-python?utm_medium=organic&utm_
    # source=google_rich_qa&utm_campaign=google_rich_qa
    if isinstance(o, np.int64):
        return int(o)
    raise TypeError


def _tup_and_byte(obj):
    """ this is used in loading """

    # convert all strings to bytes
    if isinstance(obj, (bytes)):
        return obj.decode()  # encode('utf-8')
        #return obj.encode('utf-8')

    # if this is a list of values, return list of byteified values
    if isinstance(obj, list):
        return [_tup_and_byte(item) for item in obj]

    # if this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(obj, dict):
        if "__tuple__" in obj:
            return tuple(_tup_and_byte(item) for item in obj["items"])
        else:
            return {
                _tup_and_byte(key): _tup_and_byte(val) for
                key, val in obj.items()
                }

    # if it's anything else, return it in its original form
    return obj


## FIXME: This will need to be modified to fit the new program
def _save_json(data, quiet=False):
    """
    Save PIED as json
    ## data as dict
    """
    # store params without the reference to Assembly object in params
    paramsdict = {i: j for (i, j) in data.paramsdict.items() if i != "_data"}

    # store all other dicts
    datadict = OrderedDict([\
        ("name", data.name),\
        ("__version__", data._version),\
        ("_sep", data._sep),\
        ("paramsdict", paramsdict),\
        ("_hackersonly", data._hackersonly),\
    ])

    ## FIXME: Saved as an example of how to store a more complicated data
    ##        structure
    ## save taxa
    #taxadict = OrderedDict([])
    #for key, taxon in data.taxa.items():
    #    taxadict[key] = taxon._to_fulldict()

    ## json format it using cumstom Encoder class
    fulldumps = json.dumps({
        "model": datadict,
    #    "taxa": taxadict
    },
        cls=_Encoder,
        sort_keys=False, indent=4, separators=(",", ":"),
        default=_default,
    )

    ## save to file
    modelpath = os.path.join(data.paramsdict["project_dir"],\
                                data.name + ".json")
    if not os.path.exists(data.paramsdict["project_dir"]):
        os.mkdir(data.paramsdict["project_dir"])

    ## protect save from interruption
    done = 0
    if not quiet: print("  Saving Core to {}".format(modelpath))
    while not done:
        try:
            with open(modelpath, 'w') as jout:
                jout.write(fulldumps)
            done = 1
        except (KeyboardInterrupt, SystemExit):
            print('.')
            continue


#############################
## Model Parameter Info Dicts
#############################
## FIXME: Add real parameters here
PARAMS = {
    "simulation_name" : "The name of this simulation scenario",\
    "project_dir" : "Where to save files",\
    "birth_rate" : "Speciation rate",\
    "stop_criterion" : "Whether to stop on ntaxa or time",\
    "ntaxa" : "Number of taxa to simulate if stop is `ntaxa`",\
    "time" : "Amount of time to simulate if stop is `time`",\
    "process" : "Whether to evolve `abundance` or growth `rate` via BM",\
    "ClaDS" : "Whether to allow speciation rates to change along the branches a la ClaDS",\
    "abundance_mean" : "Ancestral abundance at time 0",\
    "abundance_sigma" : "Rate at which abundance changes if process is `abundance`",\
    "growth_rate_mean" : "Ancestral population growth rate at time 0.",\
    "growth_rate_sigma" : "Rate at which growth rate changes if process is `rate`",\
    "lambda_mean" : "Ancestral speciation rate at time 0.",\
    "lambda_sigma" : "Rate at which speciation rate changes if ClaDS is True.",\
    "alpha" : "Rate shift if ClaDS is True",\
    "mutation_rate" : "Mutation rate per base per generation",\
    "sample_size" : "Number of samples to draw for calculating genetic diversity"
}


#############################
## Global error messages
#############################
BAD_PIED_NAME = """\
    No spaces or special characters of any kind are allowed in the simulation
    name. Special characters include all punctuation except dash '-' and
    underscore '_'. A good practice is to replace spaces with underscores '_'.
    An example of a good simulation name is: hawaiian_arthropods

    Here's what you put:
    {}
    """

PARAMS_EXISTS = """
    Error: Params file already exists: {}
    Use force argument to overwrite.
    """

NO_OUTDIR = """
    Error: Attempting to write params to a directory that doesn't exist - {}
    """

REQUIRE_NAME = """\
    Simulation scenario name _must_ be set. This is the first parameter in the
    params.txt file, and will be used as a prefix for output files. It should be a
    short string with no special characters, i.e., not a path (no \"/\" characters).
    If you need a suggestion, name it after the location you're working on.
    """

if __name__ == "__main__":
    pass
