#Normal Python modules
import itertools
from copy import deepcopy,copy
import numpy as np
import multiprocessing as mp
import pickle

##Fwdpy and libsequence
import fwdpy11 as fp11 #fwdpy11 
from fwdpy11.model_params import SlocusParams
import fwdpy11.wright_fisher as wf
import fwdpy11.sampling as fps
import libsequence.polytable as polyt
from libsequence.summstats import PolySIM

#My stuff
import clonal  #This is my fwdpy11 extension module for clonal evolution
import RecordStats # My recorder class

def main(rep):
    #Basic values
    rng0=fp11.GSLrng(rep+1)
    rng1=fp11.GSLrng(rep+2)
    rng2=fp11.GSLrng(rep+3)
    rng3=fp11.GSLrng(rep+4)

    #Basic values
    N_anc=4000 # Ancestral size
    N_pre=200 # Present size
    Gen_burnin = 8*N_anc # Burnin generations
    Gen_bottle = 100 #Bottle neck generations
    theta=4
    rho=1
    #rep=1
    #Initialize recorder -- every 100 generations at first
    recorder_anc = RecordStats.RecordStats(100,rng0)

    #Parameter dictionarys
    # Ancestral burnin
    p_anc = {'demography':np.array([N_anc]*Gen_burnin,dtype=np.uint32),
       'nregions':[fp11.Region(0,1,1)],
        'recregions':[fp11.Region(0,1,1)],
        'sregions':[fp11.GammaS(0,1,1,-0.05,0.3,1)],
        'rates':((11./12.)*theta/float(4*N_anc),(1./12.)*theta/float(4*N_anc),rho/float(4*N_anc))
    }

    # Bottle neck parameters
    p_bottle ={'demography':np.uint32(np.linspace(N_anc,N_pre,Gen_bottle)),
       'nregions':[fp11.Region(0,1,1)],
        'recregions':[fp11.Region(0,1,1)],
        'sregions':[fp11.GammaS(0,1,1,-0.05,0.3,1)],
        'rates':((11./12.)*theta/float(4*N_anc),(1./12.)*theta/float(4*N_anc),rho/float(4*N_anc))
    }

    #Continue constant pop for same duration as bottleneck
    p_const ={'demography':np.array([N_anc]*Gen_bottle,dtype=np.uint32),
       'nregions':[fp11.Region(0,1,1)],
        'recregions':[fp11.Region(0,1,1)],
        'sregions':[fp11.GammaS(0,1,1,-0.05,0.3,1)],
        'rates':((11./12.)*theta/float(4*N_anc),(1./12.)*theta/float(4*N_anc),rho/float(4*N_anc))
    }

    params_anc=SlocusParams(**p_anc)
    params_bottle=SlocusParams(**p_bottle)
    params_const=SlocusParams(**p_const)

    pop_anc = fp11.SlocusPop(N_anc)

    print("Burnin")

    fp11.wright_fisher.evolve(rng0,pop_anc,params_anc,recorder_anc)
    
    recorder_const_outcross = RecordStats.RecordStats(1,rng0)    

    pop_bottle_outcross = deepcopy(pop_anc) # bottle neck + outcrossing
    recorder_bottle_outcross = RecordStats.RecordStats(1,rng1) #copy(recorder_anc)

    pop_const_clonal = deepcopy(pop_anc) # constant pop size + clonal propagation
    recorder_const_clonal = RecordStats.RecordStats(1,rng2)

    pop_bottle_clonal = deepcopy(pop_anc) # bottle neck + clonal
    recorder_bottle_clonal = RecordStats.RecordStats(1,rng3)

    print("Finishing")
    fp11.wright_fisher.evolve(rng0,pop_anc,
                              params_const,recorder_const_outcross) # constant pop  + outcrossing

    fp11.wright_fisher.evolve(rng1,pop_bottle_outcross,
                              params_bottle,recorder_bottle_outcross) # bottleneck  + outcrossing

    clonal.evolve(rng2,pop_const_clonal,
                  params_const,recorder_const_clonal) # constant pop  + clonal

    clonal.evolve(rng3,pop_bottle_clonal,
                  params_bottle,recorder_bottle_clonal)# constant pop  + clonal
    print("Writing")
    recorder_dict = {"ancestral" : recorder_anc,
                    "const_outcross" : recorder_const_outcross,
                    "bottle_outcross" : recorder_bottle_outcross,
                    "const_clonal" : recorder_const_clonal,
                    "bottle_clonal" : recorder_bottle_clonal}
    
    output = open('data/test/small_clonal_grape.'+str(rep)+'.pkl', 'wb')
    pickle.dump(recorder_dict,output)
    output.close()
if __name__ == "__main__":
    #main()
    repid = [i for i in range(500)]
    P=mp.Pool()
    res=P.imap_unordered(main,repid)
    P.close()
    P.join()
