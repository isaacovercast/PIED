import datetime
import json
import logging
import momi
import numpy as np
import os
import string
import time
import tempfile

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
            ## Cast params to correct types
            if param == "project_dir":
                ## If it already exists then just inform the user that we'll be adding
                ## more simulations to the current project directory
                if " " in newvalue:
                    raise PIEDError("`project_dir` may not contain spaces. You put:\n{}".format(newvalue))
                self.paramsdict[param] = os.path.realpath(os.path.expanduser(newvalue))

                if not os.path.exists(self.paramsdict["project_dir"]):
                    os.mkdir(self.paramsdict["project_dir"])
            
            elif param in ["N_e", "tau", "epsilon"]:
                tup = tuplecheck(newvalue, dtype=int)
                if isinstance(tup, tuple):
                    self.paramsdict[param] = tup
                    if tup[0] <= 0:
                        raise PIEDError("{} values must be strictly > 0. You put {}".format(param, tup))
                else:
                    self.paramsdict[param] = tup
                    if tup <= 0:
                        raise PIEDError("{} values must be strictly > 0. You put {}".format(param, tup))
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
        npops = self.paramsdict["npops"]
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
        npops = self.paramsdict["npops"]
    
        msfs_list = []

        printstr = " Performing Simulations    | {} |"
        start = time.time()
        for i in range(nsims):
            try:
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                if not quiet: progressbar(nsims, i, printstr.format(elapsed))

                sfs_list = []
                for tidx, tau_pops in enumerate(pops_per_tau):
                    for pidx in range(tau_pops):
                        name = "pop{}-{}".format(tidx, pidx)
                        sfs_list.append(self._simulate(name,
                                                N_e=N_es[tidx],
                                                tau=taus[tidx],
                                                epsilon=epsilons[tidx]))
                msfs = multiSFS(sfs_list,\
                                sort=self._hackersonly["sorted_sfs"],\
                                proportions=self._hackersonly["proportional_msfs"])

                ## In the pipe_master model the first tau in the list is the co-expansion time
                ## If/when you get around to doing the msbayes model of multiple coexpansion
                ## pulses, then this will have to change 
                msfs.set_params(pd.Series([zeta, zeta_e, psi, taus[0], pops_per_tau, taus, epsilons, N_es],\
                                        index=["zeta", "zeta_e", "psi", "t_s", "pops_per_tau", "taus", "epsilons", "N_es"]))
                msfs_list.append(msfs)

            except KeyboardInterrupt as inst:
                print("\n    Cancelling remaining simulations")
                break
            except Exception as inst:
                LOGGER.debug("Simulation failed: {}".format(inst))
                raise PIEDError("Failed inside serial_simulate: {}".format(inst))

        if not quiet: progressbar(100, 100, " Finished {} simulations in   {}\n".format(i+1, elapsed))

        return msfs_list


    def _simulate(self, name, N_e=1e6, tau=20000, epsilon=10, verbose=False):
        model = momi.Core(N_e=N_e)
        model.add_leaf(name)
        if epsilon > 0:
            # epsilon positive, expansion
            model.set_size(name, t=tau, N=N_e/epsilon)
        elif epsilon < 0:
            # epsilon negative, bottleneck
            model.set_size(name, t=tau, N=N_e/epsilon)
        else:
            # epsilon == 0 no size change
            pass
        sampled_n_dict={name:self.paramsdict["nsamps"]}
        if verbose: print(sampled_n_dict)
        ac = model.simulate_data(length=self.paramsdict["length"],
                                num_replicates=self.paramsdict["num_replicates"],
                                recoms_per_gen=self.paramsdict["recoms_per_gen"],
                                muts_per_gen=self._sample_mu(),
                                sampled_n_dict=sampled_n_dict)
        try:
            sfs = ac.extract_sfs(n_blocks=1)
        except ValueError:
            ## If _sample_mu() returns zero, or a very small value with respect to
            ## sequence length, Ne, and tau, then you can get a case where there
            ## are no snps in the data, and constructing the sfs freaks.
            raise PIEDError("Can't extract SFS from a simulation with no variation. Check that muts_per_gen looks reasonable.")

        return sfs

    
    def simulate(self, nsims=1, ipyclient=None, quiet=False, verbose=False, force=False):
        """
        Do the heavy lifting here. 

        :param int nsims: The number of PIED codemographic simulations to
            perform.
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
        param_df = pd.DataFrame([], columns=["zeta", "psi", "pops_per_tau", "taus", "epsilons"])

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
        try:
            dat = pd.read_csv(simfile, sep=self._sep)
        except FileNotFoundError:
            dat = pd.DataFrame()

        ## sort=False suppresses a warning about non-concatenation index if
        ## SIMOUT is empty
        msfs_df = pd.DataFrame(pd.concat([x.to_dataframe() for x in msfs_list], sort=False)).fillna(0)
        dat = pd.concat([dat, result_df], sort=False)
        dat.to_csv(simfile, header=True, index=False, sep=self._sep, float_format='%.3f')


def serial_simulate(model, nsims=1, quiet=False, verbose=False):
    import os
    LOGGER.debug("Entering sim - {} on pid {}\n{}".format(model, os.getpid(), model.paramsdict))
    res = model.serial_simulate(nsims, quiet=quiet, verbose=verbose)
    LOGGER.debug("Leaving sim - {} on pid {}".format(model, os.getpid()))
    return res


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
