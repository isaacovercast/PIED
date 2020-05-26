""" the main CLI for calling PIED """

import pkg_resources
import argparse
import logging
import atexit
import time
import sys
import os

import PIED
from PIED.util import *
from PIED.parallel import *

LOGGER = logging.getLogger(__name__)


def parse_params(args):
    """ Parse the params file and return params as a dict"""

    ## check that params.txt file is correctly formatted.
    try:
        with open(args.params) as paramsin:
            plines = paramsin.read()
    except IOError as _:
        sys.exit("  No params file found")

    ## The [1:] thing just drops the first result from split which is '' anyway.
    plines = plines.split("------- ")[1:]

    ## Get params and make into a dict, ignore all blank lines
    items = [i.split("##")[0].strip() for i in plines[0].split("\n")[1:] if not i.strip() == ""]
    keys = list(PIED.Core('null', quiet=True).paramsdict.keys())
    params = {str(i):j for i, j in zip(keys, items)}

    LOGGER.debug("Got params - {}".format(params))

    return params


def getmodel(args, params):
    """ 
    loads simulation or creates a new one and set its params as
    read in from the params file. Does not launch ipcluster. 
    """

    project_dir = PIED.util._expander(params['project_dir'])
    LOGGER.debug("project_dir: {}".format(project_dir))
    sim_name = params['simulation_name']

    ## make sure the working directory exists.
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)

    data = PIED.Core(sim_name, quiet=args.quiet, verbose=args.verbose)

    ## Populate the parameters
    for param in params:
        try:
            data = set_params(data, param, params[param], quiet=args.quiet)
        except Exception as inst:
            print("Error in __main__.getmodel(): {}".format(inst))
            sys.exit(-1)

    return data


def parse_command_line():
    """ Parse CLI args."""

    ## create the parser
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 
    PIED -n data                       ## create new file called params-data.txt 
    PIED -p params-data.txt            ## run PIED with settings in params file
    """)

    ## add arguments 
    parser.add_argument('-n', metavar='new', dest="new", type=str, 
        help="create new file 'params-{new}.txt' in current directory")

    parser.add_argument('-p', metavar='params', dest="params",
        type=str,
        help="path to params file simulations: params-{name}.txt")

    parser.add_argument("-s", metavar="sims", dest="sims",
        type=int, default=0,
        help="Generate specified number of simulations")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=-1,
        help="number of CPU cores to use (Default=0=All)")

    parser.add_argument('-f', action='store_true', dest="force",
        help="force overwrite of existing data")

    parser.add_argument('-v', action='store_true', dest="verbose",
        help="do not print to stderror or stdout.")

    parser.add_argument('-q', action='store_true', dest="quiet",
        help="do not print anything ever.")

    parser.add_argument('-d', action='store_true', dest="debug",
        help="print lots more info to PIED_log.txt.")

    parser.add_argument('-V', action='version', 
        version=str(pkg_resources.get_distribution('PIED')),
        help=argparse.SUPPRESS)

    parser.add_argument("--ipcluster", metavar="ipcluster", dest="ipcluster",
        type=str, nargs="?", const="default",
        help="connect to ipcluster profile")

    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ## parse args
    args = parser.parse_args()

    if not any(x in ["params", "new"] for x in list(vars(args).keys())):
        print("\nBad arguments: must include at least one of"\
                +" `-p` or `-n`\n")
 
        #parser.print_help()
        sys.exit(1)

    return args


def do_sims(data, args):
    ## if ipyclient is running (and matched profile) then use that one
    if args.ipcluster:
        ipyclient = ipp.Client(cluster_id=args.ipcluster)
        data._ipcluster["cores"] = len(ipyclient.ids)
        if not args.quiet:
            print("    Attached to cluster {} w/ {} engines.".format(args.ipcluster, data._ipcluster["cores"]))

    ## if not then we need to register and launch an ipcluster instance
    elif args.cores >= 0:
        ## set CLI ipcluster terms
        ipyclient = None
        data._ipcluster["cores"] = args.cores if args.cores else detect_cpus()
        data._ipcluster["engines"] = "Local"
        ## register to have a cluster-id with "ip- name"
        data = register_ipcluster(data)
        ipyclient = get_client(**data._ipcluster)
        print(cluster_info(ipyclient))
    ## If args.cores is negative then don't use ipcluster
    else:
        ipyclient = None
        if not args.quiet: print("    Parallelization disabled.")
    try:
        ## Do stuff here
        data.simulate(
            nsims=args.sims,
            ipyclient=ipyclient,
            quiet=args.quiet,
            verbose=args.verbose,
            force=args.force)
    except KeyboardInterrupt as inst:
        print("\n  Keyboard Interrupt by user")
        LOGGER.info("assembly interrupted by user.")
    except PIEDError as inst:
        LOGGER.error("PIEDError: %s", inst)
        print("\n  Encountered an error (see details in ./PIED_log.txt)"+\
              "\n  Error summary is below -------------------------------"+\
              "\n{}".format(inst))
    except Exception as inst:
        LOGGER.error(inst)
        print("\n  Encountered an unexpected error (see ./PIED_log.txt)"+\
              "\n  Error message is below -------------------------------"+\
              "\n{}".format(inst))
    finally:
        try:
            ## can't close client if it was never open
            if ipyclient:
                if args.verbose: print("Clean up ipcluster {}".format(ipyclient))
                ## send SIGINT (2) to all engines
                try:
                    ipyclient.abort(block=False)
                    time.sleep(1)
                    for engine_id, pid in list(data._ipcluster["pids"].items()):
                        LOGGER.debug("  Cleaning up ipcluster engine/pid {}/{}".format(engine_id, pid))
                        if ipyclient.queue_status()[engine_id]["tasks"]:
                            os.kill(pid, 2)
                            LOGGER.info('interrupted engine {} w/ SIGINT to {}'\
                                    .format(engine_id, pid))
                    time.sleep(1)
                except ipp.NoEnginesRegistered as inst:
                    LOGGER.debug("No engines registered: {}".format(inst))

                ## if CLI, stop jobs and shutdown. Don't use _cli here 
                ## because you can have a CLI object but use the --ipcluster
                ## flag, in which case we don't want to kill ipcluster.
                if 'PIED-cli' in data._ipcluster["cluster_id"]:
                    LOGGER.info("  shutting down engines")
                    ipyclient.shutdown(hub=True, block=False)
                    ipyclient.close()
                    LOGGER.info("  finished shutdown")
                else:
                    if not ipyclient.outstanding:
                        ipyclient.purge_everything()
                    else:
                        ## nanny: kill everything, something bad happened
                        ipyclient.shutdown(hub=True, block=False)
                        ipyclient.close()
                        print("\nwarning: ipcluster shutdown and must be restarted")

        ## if exception is close and save, print and ignore
        except Exception as inst2:
            print("warning: error during shutdown:\n{}".format(inst2))
            LOGGER.error("shutdown warning: %s", inst2)


def main():
    """ main function """
    PIED.__interactive__ = 0  ## Turn off API output

    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    if not args.quiet: print(PIED_HEADER)

    ## Turn the debug output written to PIED_log.txt up to 11!
    ## Clean up the old one first, it's cleaner to do this here than
    ## at the end (exceptions, etc)
    if os.path.exists(PIED.__debugflag__):
        os.remove(PIED.__debugflag__)

    if args.debug:
        if not args.quiet: print("\n  ** Enabling debug mode **\n")
        PIED._debug_on()
        atexit.register(PIED._debug_off)

    ## create new paramsfile if -n
    if args.new:
        ## Create a tmp assembly, call write_params to make default params.txt
        try:
            tmpassembly = PIED.Core(args.new, quiet=True)
            tmpassembly.write_params("params-{}.txt".format(args.new), outdir="./", 
                                     force=args.force)
        except Exception as inst:
            print("Error creating new params file: {}".format(inst))
            sys.exit(2)

        print("\n  New file 'params-{}.txt' created in {}\n".\
               format(args.new, os.path.realpath(os.path.curdir)))
        sys.exit(2)


    ## if params then must provide action argument with it
    if args.params:
        if not any([args.sims]):
            print(PIED_USAGE)
            sys.exit(2)

    if not args.params:
        if any([args.sims]):
            print(PIED_USAGE)
            sys.exit(2)

    ## Log the current version. End run around the LOGGER
    ## so it'll always print regardless of log level.
    with open(PIED.__debugfile__, 'a') as logfile:
        logfile.write(PIED_HEADER)
        logfile.write("\n  Begin run: {}".format(time.strftime("%Y-%m-%d %H:%M")))
        logfile.write("\n  Using args {}".format(vars(args)))
        logfile.write("\n  Platform info: {}\n".format(os.uname()))

    ## create new simulation model
    if args.params:
        params = parse_params(args)
        LOGGER.debug("params - {}".format(params))

        ## launch or load simulationwith custom profile/pid
        data = getmodel(args, params)

        ## Generate numerous simulations
        if args.sims:
            ## Only blank the log file if we're actually going to do real
            ## work. In practice the log file should never get this big.
            if os.path.exists(PIED.__debugfile__):
                if os.path.getsize(PIED.__debugfile__) > 50000000:
                    with open(PIED.__debugfile__, 'w') as clear:
                        clear.write("file reset")

            #if not args.quiet:
            #    print("\n    {}".format(data))

            try:
                do_sims(data, args)
            except Exception as inst:
                print("  Unexpected error - {}".format(inst))


BAD_PARAMETER_ERROR = """
    Malformed params file: {}
    Bad parameter {} - {}
    {}"""

PIED_HEADER = \
"\n -------------------------------------------------------------"+\
"\n  PIED [v.{}]".format(PIED.__version__)+\
"\n  PhIgurE out an acronymD"+\
"\n -------------------------------------------------------------"


PIED_USAGE = """
    Must provide action argument along with -p argument for params file. 
    e.g., PIED -p params-test.txt -s 10000           ## run 10000 simulations
    """

if __name__ == "__main__": 
    main()
