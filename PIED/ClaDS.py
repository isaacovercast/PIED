import numpy as np
import pandas as pd
import toytree

from toytree.utils import TreeError

## Only implements the 'turnover' method for mu, so no death parameter
def ClaDS(ntips=10,
            time=4,
            b=1,
            stop="taxa",
            seed=None,
            retain_extinct=False,
            ClaDS_sigma=0.2,
            ClaDS_alpha=0.7,
            mu_0=0.5,
            verbose=False):

    if stop not in ["taxa", "time"]:
        raise TreeError("stop must be either 'taxa' or 'time'")

    if seed:
        np.random.seed(seed)

    tre = toytree.tree()
    rt = tre.treenode.add_child(name="0", dist=0)
    rt.add_feature("tdiv", 0)
    rt.add_feature("lambda_", b)
    rt.add_feature("mu", b*mu_0)

    taxa_stop = ntips
    time_stop = time

    # Counters for extinctions, total events, and time
    ext = 0
    evnts = 0
    t = 0

    while(1):
        tips = tre.treenode.get_leaves()

        lambs = [tip.lambda_ for tip in tips]
        mus = [tip.mu for tip in tips]
        # Run a horse race for all lineages, smallest time sampled wins
        # This is exactly equal to the way ClaDS does it, this way makes
        # more sense to me. See jupyter-notebooks/misc_util.ipynb.
        ts = np.random.exponential(1/(np.array(lambs+mus)))
        idx = np.where(ts == ts.min())[0][0]

        death = False
        # Figure out if it's a birth or death event
        # Shift the idx down by the number of tips and set the flag
        if idx > len(lambs)-1:
            idx = idx - len(lambs)
            death = True

        dt = ts.min()
        sp = tips[idx]

        t = t + dt
        evnts += 1


        if not death:
            # birth event
            c1 = sp.add_child(name=str(t)+"-1", dist=0)
            c2 = sp.add_child(name=str(t)+"-2", dist=0)
            for c in [c1, c2]:
                c.add_feature("tdiv", t)
                c.add_feature("lambda_", np.random.lognormal(np.log(sp.lambda_ * ClaDS_alpha), ClaDS_sigma))
                c.add_feature("mu", c.lambda_ * mu_0)
        else:
            # extinction event
            if retain_extinct:
                # Deletes the tip node, but retains the node leading
                # to the sister
                sp.delete(preserve_branch_length=True,\
                            prevent_nondicotomic=False)
                tre = _prune(tre)
            else:
                try:
                    tips.remove(sp)
                    tre.treenode.prune(tips, preserve_branch_length=True)
                except TreeError:
                    # The root lineage can't go extinct
                    pass
            ext += 1

        # Update branch lengths
        tips = tre.treenode.get_leaves()
        for x in tips:
            x.dist += dt

        # Check stopping criterion
        tips = tre.treenode.get_leaves()
        done = False
        if stop == "taxa":
            if len(tips) >= taxa_stop:
                done = True
        elif stop == "time":
            if t >= time_stop:
                done = True
        elif len(tips) == 0:
            print("All lineages extinct")
            done = True
        if done:
            if verbose:
                print("Birth events {}".format(evnts))
                print("Extinctions {} (per birth {})".format(ext, ext/evnts))
            for i, t in enumerate(tips[::-1]):
                t.name = "r{}".format(i)
            tre._coords.update()
            return tre


def _prune(tre, verbose=False):
    "Helper function for recursively pruning extinct branches in bd trees"
    ttree = tre.copy()
    tips = ttree.treenode.get_leaves()

    if np.any(np.array([x.height for x in tips]) > 0):
        for t in tips:
            if not np.isclose(t.height, 0):
                if verbose: print("Removing node/height {}/{}".format(t.name,\
                                                                      t.height))
                t.delete(prevent_nondicotomic=False)
                ttree = _prune(ttree)
    return ttree
