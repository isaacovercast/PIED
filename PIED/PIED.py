import ast
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
from scipy.stats import hmean

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
                       ("ClaDS_sigma", 0.1),
                       ("ClaDS_alpha", 0.1),
                       ("sequence_length", 500),
                       ("mutation_rate", 1e-5),
                       ("sample_size", 10),
                       ("abundance_scaling", "None")
        ])

        ## Separator to use for reading/writing files
        self._sep = " "

        ## A dictionary for storing taxon specific information. This dictionary
        ## is populated when empirical data is loaded. If no per taxon info
        ## is present then we sample from the priors as specified in the params
        ## file.
        self.taxa = {}

        ## elite hackers only internal dictionary, normally you shouldn't mess with this
        ## max_theta - The max value of theta, before we consider the sim to be
        ##             be degenerate and bail out.
        ## max_t -  The maximum amount of elapsed time in a given simulation
        ##          before we consider it degenerate and bail out. In units mya.
        ##          This is a somewhat unsatisfying solution to a problem that can
        ##          arise in the clads model. Occasionally as the speciation rate
        ##          becomes lower and lower, the waiting time becomes higher and
        ##          higher (astronomical), and eventually it's so long that
        ##          abundances become infinite and it flips the board.
        ## harmonic_mean - Whether to take the harmonic mean of changing abundance
        ##                 through time, or simply use current abundance as Ne
        self._hackersonly = dict([
                        ("max_theta", 0.7),
                        ("max_t", 3000),
                        ("harmonic_mean", True),
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
                self.paramsdict[param] = newvalue
            elif param == "process":
                if newvalue not in ["abundance", "rate"]:
                    raise PIEDError("Bad parameter: `process` must be "\
                                    + "`abundance` or `rate`.")
                self.paramsdict[param] = newvalue
            elif param == "ClaDS":
                try:
                    self.paramsdict[param] = ast.literal_eval(newvalue)
                except ValueError:
                    raise PIEDError("Bad parameter: `ClaDS` must be `True` or "\
                                    + "`False`.")
            elif param == "abundance_scaling":
                if newvalue in ["log", "ln", "None"]:
                    self.paramsdict[param] = newvalue
                else:
                    try:
                        self.paramsdict[param] = float(newvalue)
                    except ValueError:
                        raise PIEDError("Bad parameter: `abundance_scaling` must be "\
                                        +"'None', 'log', 'ln' or a ratio as a float value.")
            ## All strictly positive float parameters
            elif param in ["abundance_sigma", "ClaDS_sigma", "ClaDS_alpha",
                            "birth_rate", "time", "growth_rate_sigma",
                            "mutation_rate"]:
                tup = tuplecheck(newvalue, dtype=float)
                msg = "Bad parameter: `{}` must be strictly positive.".format(param)
                if isinstance(tup, tuple):
                    if tup[0] <= 0:
                        raise PIEDError(msg)
                elif tup <= 0:
                        raise PIEDError(msg)
                self.paramsdict[param] = tup
            elif param in ["ntaxa", "abundance_mean", "sequence_length",
                            "sample_size"]:
                tup = tuplecheck(newvalue, dtype=int)
                msg = "Bad parameter: `{}` must be strictly positive.".format(param)
                if isinstance(tup, tuple):
                    if tup[0] <= 0:
                        raise PIEDError(msg)
                    elif tup <= 0:
                        raise PIEDError(msg)
                self.paramsdict[param] = tup
            ## Growth rate mean can be zero, no problem
            elif param == "growth_rate_mean":
                tup = tuplecheck(newvalue, dtype=float)
                self.paramsdict[param] = tup
            else:
                self.paramsdict[param] = newvalue
        except Exception as inst:
            LOGGER.debug("Error setting parameter: {} {}".format(param, newvalue))
            raise


    ## Getting parameters header and parameters carves off
    ## the simulation name and the project directory
    ## Slice off the first two params which are name and project_dir
    def _get_params_header(self):
        return list(self.paramsdict.keys())[2:]


    def _get_params_values(self):
        ## This is only ever called when writing out to the simout file
        ## so make str here to handle int/float params
        return [str(x) for x in list(self.paramsdict.values())[2:]]


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


    def get_params(self, verbose=False):
        """
        A convenience function for getting nicely formatted params in API mode.

        :return: A string of all the params ready to be printed.
        """
        tf = tempfile.NamedTemporaryFile()
        self.write_params(outfile=tf.name, force=True)
        dat = open(tf.name).read()
        if verbose: print(dat)
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
                LOGGER.error("Failed inside parallel_simulate -\n{}".format(inst))
                ## Don't let one bad apple spoil the bunch,
                ## so keep trying through the rest of the asyncs
        if len(result_list) < nsims:
            print("\n   One or more simulations failed. Check PIED_log.txt for details.\n")

        return result_list
    
    
    def serial_simulate(self, nsims=1, quiet=False, verbose=False):
        import pandas as pd

        tree_list = []

        printstr = " Performing Simulations    | {} |"
        start = time.time()
        fail_list = []
        for i in range(nsims):
            try:
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                if not quiet: progressbar(nsims, i, printstr.format(elapsed))

                ## Returns the tree and the formatted results
                _, res = self._simulate()
                tree_list.append(res)

            except KeyboardInterrupt as inst:
                print("\n    Cancelling remaining simulations")
                break
            except Exception as inst:
                LOGGER.error("\n\tSimulation failed: {}".format(inst))
                fail_list.append(inst)

        if not quiet: progressbar(100, 100, " Finished {} simulations in   {}\n".format(i+1, elapsed))
        if len(fail_list) > 0:
            print("\n  {} failed simulations inside serial_simulate. See PIED_log.txt for details.\n".format(len(fail_list)))

        return tree_list


    def _simulate(self, verbose=False):
        
        # Each species has an attribute for each feature in this dictionary.
        #  abundance - Abundance of the species
        #  r - Rate at which abundance changes, can be negative
        #  trait - This is a random trait value that evolves by BM
        #  lambda - Per lineage speciation rate
        feature_dict = {"abundance":{"sigma":self.paramsdict["abundance_sigma"], "zbar_0":self.paramsdict["abundance_mean"], "log":True, "dtype":"int"},
                          "r":{"sigma":self.paramsdict["growth_rate_sigma"], "zbar_0":self.paramsdict["growth_rate_mean"], "log":False, "dtype":"float"},
                          "trait":{"sigma":2, "zbar_0":0, "log":False, "dtype":"float"},
                         }

        tre = toytree.tree()
        ## Items in the feature_dict all evolve via BM in most cases
        for fname, fdict in feature_dict.items():
            tre.treenode.add_feature(fname, fdict["zbar_0"])
        ## Add the default lambda_ which is only used in ClaDS
        tre.treenode.add_feature("lambda_", self.paramsdict["birth_rate"])
        ## Add a list to track abundances through time
        tre.treenode.add_feature("abunds", [])

        taxa_stop = self.paramsdict["ntaxa"]
        time_stop = self.paramsdict["time"]

        ext = 0
        evnts = 0
        t = 0
        while(1):

            ## Get list of extant tips
            tips = tre.treenode.get_leaves()

            # Sample time interval
            if self.paramsdict["ClaDS"]:
                lambs = np.array([tip.lambda_ for tip in tips])
                # Run a horse race for all lineages, smallest time sampled wins
                # This is exactly equal to the way ClaDS does it, this way makes
                # more sense to me. See jupyter-notebooks/misc_util.ipynb.
                ts = np.random.exponential(1/lambs)
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
                ## If abundance == 1 then the call to randint will raise.
                ## This _used_ to happen, but then I think I fixed it. You can
                ## remove this the next time you feel like it.
                if sp.abundance == 1:
                    raise PIEDError("Sp w/ abundance=1. This should never happen.")
                else:
                    raise PIEDError("Abundance exceeds int64: {}".format(sp.abundance))

            c1.add_feature("abundance", abund)
            c2.add_feature("abundance", sp.abundance-abund)

            ## Set new species features
            for c in [c1, c2]:
                c.abunds = [c.abundance]
                c.r = sp.r
                c.trait = sp.trait
                if self.paramsdict["ClaDS"]:
                    c.lambda_ = np.random.lognormal(np.log(sp.lambda_\
                                                     * self.paramsdict["ClaDS_alpha"]),\
                                                       self.paramsdict["ClaDS_sigma"])
                else:
                    c.lambda_ = sp.lambda_

            # Update branch length, evolve features, and update abundance
            tips = tre.treenode.get_leaves()
            for x in tips:
                x.dist += dt
                for fname, fdict in feature_dict.items():
                    # Update 'r' and 'trait', but skip 'abundance'
                    # if doing the growth rate process.
                    if self.paramsdict["process"] == "rate" and\
                        fname == "abundance": continue
                    else:
                        x.add_feature(fname, _bm(mean=getattr(x, fname),
                                                sigma=fdict["sigma"],
                                                dt=dt,
                                                log=fdict["log"],
                                                dtype=fdict["dtype"]))
                if self.paramsdict["process"] == "rate":
                    # Apply the population size change
                    try:
                        x.abundance = int(x.abundance * (np.exp(x.r*dt)))
                    except Exception as inst:
                        raise PIEDError("Abundance is too big: {}".format(x.abundance))
                ## Record abundance through time
                x.abunds.append(x.abundance)

            # Check boundary conditions for abundance and lambda
            extinct = []
            for x in tips:
                # Check extinction and prune out extinct and any resulting hanging branches
                # If you're the only individual left, you have nobody to mate with so you die.
                if x.abundance <= 1:
                    extinct.append(x)
                    ext += 1

            tips = [x for x in tips if x not in extinct]
            if len(tips) == 0:
                # Can't prune to an empty tree, so simply return a new empty tree
                tre = toytree.tree()
            else:
                tre.treenode.prune(tips, preserve_branch_length=True)

            ## This shouldn't be necessary now that we're using the treenode.prune() method
            #tre = _prune(tre)

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
                raise PIEDError("All lineages extinct")
            elif t >= self._hackersonly["max_t"]:
                raise PIEDError("Max time exceded")
            if done:
                if verbose:
                    print("ntips {}".format(len(tips)))
                    print("time {}".format(t))
                    print("Birth events {}".format(evnts))
                    print("Extinctions (per birth) {} ({})".format(ext, ext/evnts))
                tre._coords.update()

                ## Scale abundance to Ne for each tip. The call to _scale_abundance
                ## returns the tip with Ne/Nes set
                tips = [_scale_abundance(t, self.paramsdict["abundance_scaling"]) for t in tips]

                ## Test for 'reasonable' theta values. If theta gets too large
                ## two things happen: biologically meaningless pi values, and
                ## runtime of msprime _explodes_.
                thetas = np.array([x.Ne * self.paramsdict["mutation_rate"] for x in tips])
                if np.any(thetas > self._hackersonly["max_theta"]):
                    raise PIEDError(MAX_THETA_ERROR.format(self._hackersonly["max_theta"],\
                                                            np.max(thetas)))

                ## Relabel tips to have reasonable names
                for i, tip in enumerate(tips[::-1]):
                    tip.name = "r{}".format(i)
                    ## Scale abundance through time to Ne by some fashion
                    tip.pi = nucleotide_diversity(self.paramsdict,
                                                    node=tip,
                                                    harmonic=self._hackersonly["harmonic_mean"])
                ## Build the output list to return which will include
                ##  * parameters of the model
                ##  * observed ntaxa and time, and calculated extinction rate
                ##  * All the data at the tips including abundance, pop growth rate, and speciation rate
                ##  * The newick formatted tree
                res = self._get_params_values()
                res.extend([str(len(tips)), str(t), str(ext/evnts)])
                dat = ["{}:{}:{}:{}:{}".format(x.name, x.abundance, x.pi, x.r, x.lambda_) for x in tips]
                dat = ",".join(dat)
                res.append(dat)
                res.append(tre.write())
                return tre, res

    
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
        if force:
            ## Prevent from shooting yourself in the foot with -f
            try:
                os.rename(simfile, simfile+".bak")
            except FileNotFoundError:
                ## If the simfile doesn't exist catch the error and move on
                pass

        if (not os.path.exists(simfile)) or force:
            params = self._get_params_header()
            params.extend(["obs_ntaxa", "obs_time", "turnover_rate"])
            params.append("data")
            params.append("tree")
            with open(simfile, 'w') as simout:
                simout.write(" ".join(params) + "\n")

        with open(simfile, 'a') as output:
            ## Don't write a newline if all simulations failed
            if result_list:
                output.write("\n".join([" ".join(x) for x in result_list]) + "\n")


def serial_simulate(model, nsims=1, quiet=False, verbose=False):
    import os
    LOGGER.debug("Entering sim - {} on pid {}\n{}".format(model, os.getpid(), model.paramsdict))
    res = model.serial_simulate(nsims, quiet=quiet, verbose=verbose)
    LOGGER.debug("Leaving sim - {} on pid {}".format(model, os.getpid()))
    return res


###########################
## Random utility functions
###########################
def _scale_abundance(tip, abundance_scaling):
    if abundance_scaling == "None":
        tip.Ne = tip.abundance
        tip.Nes = tip.abunds
    elif abundance_scaling == "ln":
        tip.Ne = np.log(tip.abundance)
        tip.Nes = np.log(tip.abunds)
    elif abundance_scaling == "log":
        tip.Ne = np.log10(tip.abundance)
        tip.Nes = np.log10(tip.abunds)
    else:
        tip.Ne = tip.abundance * abundance_scaling
        tip.Nes = [x * abundance_scaling for x in tip.abunds]
    return tip


def nucleotide_diversity(paramsdict, node, harmonic=False):
    if harmonic:
        Ne = hmean(node.Nes)
    else:
        Ne = node.Ne
    ts = msprime.simulate(sample_size=paramsdict["sample_size"],
                            Ne=Ne,
                            length=paramsdict["sequence_length"],
                            mutation_rate=paramsdict["mutation_rate"])
    ## By default tskit.diversity() is per base
    return ts.diversity()


def nucleotide_diversity_ILS(paramsdict, tree, debug=False):
    generation_time = 1
    population_configurations = []
    demographic_events = []
    msp_idx = 0
    ## Traverse the tree generating population configurations for tips and
    ## nodes. Kind of annoying.
    for i, node in enumerate(tree.treenode.traverse("postorder")):
        children = node.get_descendants()
        if len(children) == 0:
            pop = msprime.PopulationConfiguration(sample_size=5, initial_size=node.Ne)
            population_configurations.append(pop)
            node.add_feature("msprime_idx", msp_idx)
            if debug: print("I'm a tip {} - {}".format(node.idx, node.msprime_idx))
            msp_idx += 1
        else:
            chidx = [c.msprime_idx for c in children]
            gens = node.height * 1e6 / generation_time
            mig = msprime.MassMigration(time=gens, source=chidx[0], dest=chidx[1])
            demographic_events.append(mig)
            node.add_feature("msprime_idx", chidx[1])
            if debug: print("I'm a node {} ({}) [gens {}]- {}".format(node.idx, node.msprime_idx, gens,\
                                                            " ".join(map(lambda x: str(x), chidx))))
    ## Sort the demographic events
    demographic_events = sorted(demographic_events, key=lambda x: x.time)

    if debug:
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            demographic_events=demographic_events)
        dd.print_history()

    ## Do the simulation
    ts = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        length=paramsdict["sequence_length"],
        mutation_rate=paramsdict["mutation_rate"])

    pop_inds = {}
    for pop in ts.populations():
        pop_inds[pop.id] = ts.samples(pop.id)
    pis = ts.diversity(list(pop_inds.values()))

    for pi, leaf in zip(pis, tree.treenode.get_leaves()):
        leaf.pi = pi

    return ts, tree


## Brownian motion function
def _bm(mean, sigma, dt, log=True, dtype="int"):
    ret = 0
    mean = np.float(mean)
    if dtype in ["int", int]:
        # Avoid log(1)
        if not log or mean <= 1:
            ret = np.int(np.round(np.random.normal(mean, sigma*dt)))
        else:
            ret = np.int(np.round(np.exp(np.random.normal(np.log(mean), sigma*dt))))

    elif dtype in ["float", float]:
        ret = np.random.normal(mean, sigma*dt)
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
    "ClaDS_sigma" : "Rate at which speciation rate changes if ClaDS is True",\
    "ClaDS_alpha" : "Rate shift if ClaDS is True",\
    "sequence_length" : "Length of the genomic region simulated, in base pairs",\
    "mutation_rate" : "Mutation rate per base per generation",\
    "sample_size" : "Number of samples to draw for calculating genetic diversity",\
    "abundance_scaling" : "Scaling abundance to Ne. Can be None, log, ln or a ratio"
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

MAX_THETA_ERROR = """\
    One or more Theta values has exceeded the max theta cutoff: {}

    Your max theta: {}

    This can happen if mutation rate is too high or if Ne gets too large, which
    can be a function of abundances wandering too fast in the 'abundance'
    process, or (more likely) growth rate getting too big in the 'rate' process.
    """ 


if __name__ == "__main__":
    pass
