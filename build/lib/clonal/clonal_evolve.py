import numpy as np
import fwdpy11
from .cl_evolve import evolve_singlepop_regions_cpp

def evolve(rng,pop,params,recorder=None):
    """
    Evolve a population clonally

    :param rng: An instance of :class:`fwdpy11.fwdpy11_types.GSLrng`
    :param pop: An instance of :class:`fwdpy11.fwdpy11_types.SlocusPop`
    :param params: An instance of :class:`fwdpy11.model_params.SlocusParams`
    :param recorder: (None) A temporal sampler/data recorder.

    .. note::
        If recorder is None, then :class:`fwdpy11.temporal_samplers.RecordNothing` will be used.

    """
    import warnings
    #Test parameters while suppressing warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #Will throw exception if anything is wrong:
        params.validate()

    from fwdpy11.internal import makeMutationRegions,makeRecombinationRegions
    mm=makeMutationRegions(params.nregions,params.sregions)
    rm=makeRecombinationRegions(params.recregions)

    if recorder is None:
        from fwdpy11.temporal_samplers import RecordNothing
        recorder = RecordNothing()

    evolve_singlepop_regions_cpp(rng,pop,params.demography,
                params.mutrate_n,params.mutrate_s,params.recrate,mm,rm,
                params.gvalue,recorder,params.pself,
                params.prune_selected)
