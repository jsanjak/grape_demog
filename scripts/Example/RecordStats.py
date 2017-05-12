import numpy as np
import itertools
import fwdpy11 as fp11 #fwdpy11 
from fwdpy11.model_params import SlocusParams
import fwdpy11.wright_fisher as wf
import fwdpy11.sampling as fps
import libsequence.polytable as polyt
from libsequence.summstats import PolySIM

class RecordStats:
    """
    Record some basic stats about the population
    including average fitness
    """
    def __init__(self,ngen,rng):
        self.ngen=ngen
        self.generation = [] 
        self.N = []
        self.rng=rng
        self.relative_load = []
        self.segregating_load = []
        self.fixed_load = []
        self.fixed_deleterious = []
        self.fixed_neutral = []
        self.mean_deleterious_per_diploid = []
        self.mean_neutral_per_diploid = []
        self.cumulative_deleterious_frequency = []
        self.cumulative_neutral_frequency = []
        self.neutral_tajd = []
        self.total_tajd = []
        self.neutral_pi = []
        self.total_pi = []
        self.neutral_hprime = []
        self.total_hprime = []
    
    
    def __call__(self,pop):
        """
        The call operator  will be given with the whole population.
        One can then operate in a read-only way, no copies made.
        """
        if ( pop.generation % self.ngen == 0 ):
            #first get necessary information from the population
            N = len(pop.diploids) # Population size
            w = [ ind.w for ind in pop.diploids ] # fitness
            fixed_s = [mut.s for mut in pop.fixations if not mut.neutral]
            #Frequency of deleterious mutations
            cumulative_deleterious_frequency = np.sum([ j/float(2*N) for i,j in 
                                                       zip(pop.mutations,pop.mcounts) if 
                                                       not i.neutral and j > 0 ])
            cumulative_neutral_frequency = np.sum([ j/float(2*N) for i,j in 
                                                   zip(pop.mutations,pop.mcounts) if 
                                                   i.neutral and j > 0 ])
            #How many of each type of mutation does each diploid have?
            mean_deleterious_per_diploid = np.mean([len(pop.gametes[ind.first].smutations) + 
                                       len(pop.gametes[ind.second].smutations) for 
                                       ind in pop.diploids ])
            mean_neutral_per_diploid = np.mean([len(pop.gametes[ind.first].mutations) + 
                                       len(pop.gametes[ind.second].mutations) for 
                                       ind in pop.diploids])  
            #How many fixed mutations are there?
            fixed_deleterious = len(fixed_s)
            fixed_neutral = [mut.neutral for mut in pop.fixations].count(True)
            relative_load = 1 - np.mean(w)/np.max(w)
            segregating_load = 1 - np.mean(w)
            fixed_load = 1 - np.prod([1+2*s for s in fixed_s])

            #Statistics on population samples:
            samp = fps.sample_separate(self.rng,pop,100,False)
            neutral_sample = polyt.SimData([str2byte(mut,'utf-8') 
                                for mut in samp[0]])
            combined_sample = polyt.SimData([str2byte(mut,'utf-8') 
                                for mut in 
                                 list(itertools.chain.from_iterable(samp))])
            ps_neutral = PolySIM(neutral_sample)
            ps_combined = PolySIM(combined_sample)
            neutral_tajd = ps_neutral.tajimasd()
            total_tajd = ps_combined.tajimasd()
            neutral_pi = ps_neutral.thetapi()
            total_pi = ps_combined.thetapi()
            neutral_hprime = ps_neutral.hprime()
            total_hprime = ps_combined.hprime()            
            #Now append these to the recorder
            self.generation.append(pop.generation)
            self.N.append(N)
            self.relative_load.append(relative_load)
            self.segregating_load.append(segregating_load)
            self.fixed_load.append(fixed_load)
            self.fixed_deleterious.append(fixed_deleterious)
            self.fixed_neutral.append(fixed_neutral)
            self.mean_deleterious_per_diploid.append(mean_deleterious_per_diploid)
            self.mean_neutral_per_diploid.append(mean_neutral_per_diploid)
            self.cumulative_deleterious_frequency.append(cumulative_deleterious_frequency)
            self.cumulative_neutral_frequency.append(cumulative_neutral_frequency)
            self.neutral_tajd.append(neutral_tajd)
            self.total_tajd.append(total_tajd)
            self.neutral_pi.append(neutral_pi)
            self.total_pi.append(total_pi)
            self.neutral_hprime.append(neutral_hprime)
            self.total_hprime.append(total_hprime)
            