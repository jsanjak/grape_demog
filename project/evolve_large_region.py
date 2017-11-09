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

def makeParams(demog,theta,rho,Ne,h):
    p = {'demography':demog,
         'nregions':[fp11.Region(0,1,1)],
         'recregions':[fp11.Region(0,1,1)],
         'sregions':[fp11.GammaS(0,1,1,-0.05,0.3,h)],
         'rates':((11./12.)*theta/float(4*Ne),(1./12.)*theta/float(4*Ne),rho/float(4*Ne))} 
    return(p)

def main(rep):
    #Basic values
    rng0=fp11.GSLrng(rep+1)
    rng1=fp11.GSLrng(rep+2)
    rng2=fp11.GSLrng(rep+3)
    rng3=fp11.GSLrng(rep+4)
    rng4=fp11.GSLrng(rep+5)
    rng5=fp11.GSLrng(rep+6)
    rng6=fp11.GSLrng(rep+7)
    rng7=fp11.GSLrng(rep+8)

    #Basic values
    N_anc=4000 # Ancestral size
    N_pre=200 # Present size
    Gen_burnin = 8*N_anc # Burnin generations
    Gen_decl = 73 #Population decline pre-domestication
    Gen_domest = 27 #Domestication
    Gen_noneq = Gen_decl + Gen_domest # Total non-equilibrium time
    theta=4000 # 4 per kb, i.e. 40 = 10kb, 400=100kb, 4000=1MB
    rho=800 # .8 per kb, i.e. 8 = 10kb, 80=100kb, 800=1MB
    #rep=1
    
    for gen_model in ["additive","recessive"]:
        if (gen_model == "additive"):
            h = 1
        if (gen_model == "recessive"):
            h = 0    
        #Initialize recorder -- every 101 generations at first
        recorder_anc = RecordStats.RecordStats(100,rng0)

        #Demographies
        anc = np.array([N_anc]*Gen_burnin,dtype=np.uint32)
        lind = np.uint32(np.linspace(N_anc,N_pre,Gen_noneq))
        expd = np.uint32(np.exp(np.linspace(np.log(N_anc),np.log(N_pre),Gen_noneq)))
        discd = np.array([N_pre]*Gen_noneq,dtype=np.uint32)
        const = np.array([N_anc]*Gen_noneq,dtype=np.uint32)
        
        #Parameter dictionaries
        # Ancestral burnin
        p_anc = makeParams(anc,theta,rho,N_anc,h)
        
        #Linear decline
        p_lin_decl = makeParams(lind[0:Gen_decl],theta,rho,N_anc,h) #73 gens
        p_lin_domest = makeParams(lind[Gen_decl:Gen_noneq],theta,rho,N_anc,h) #27 gens
        
        # Exponential population decline
        p_exp_decl = makeParams(expd[0:Gen_decl],theta,rho,N_anc,h)
        p_exp_domest = makeParams(expd[Gen_decl:Gen_noneq],theta,rho,N_anc,h)

        # Discrete population reduction
        p_disc_decl = makeParams(discd[0:Gen_decl],theta,rho,N_anc,h)
        p_disc_domest = makeParams(discd[Gen_decl:Gen_noneq],theta,rho,N_anc,h)

        #Continue constant pop for same duration as population decline
        p_const_decl = makeParams(const[0:Gen_decl],theta,rho,N_anc,h)
        p_const_domest = makeParams(const[Gen_decl:Gen_noneq],theta,rho,N_anc,h)

        params_anc=SlocusParams(**p_anc)
        
        params_lin_decl=SlocusParams(**p_lin_decl)
        params_lin_domest=SlocusParams(**p_lin_domest)
        
        params_exp_decl=SlocusParams(**p_exp_decl)
        params_exp_domest=SlocusParams(**p_exp_domest)
        
        params_disc_decl=SlocusParams(**p_disc_decl)
        params_disc_domest=SlocusParams(**p_disc_domest)
        
        params_const_decl=SlocusParams(**p_const_decl)
        params_const_domest=SlocusParams(**p_const_domest)

        pop_anc = fp11.SlocusPop(N_anc)
        
        print("Burnin")

        fp11.wright_fisher.evolve(rng0,pop_anc,params_anc,recorder_anc)

        
        # Use ancestral pop as the constant pop
        pop_lin_decl = deepcopy(pop_anc) # linear + outcrossing
        pop_exp_decl = deepcopy(pop_anc) # exponential + outcrossing
        pop_disc_decl = deepcopy(pop_anc) # discrete + outcrossing
    
        recorder_const_decl = RecordStats.RecordStats(1,rng0)    
        recorder_lin_decl = RecordStats.RecordStats(1,rng1) #copy(recorder_anc)
        recorder_exp_decl = RecordStats.RecordStats(1,rng2)
        recorder_disc_decl = RecordStats.RecordStats(1,rng3)
        
        

        print("Pre-domestication")
        fp11.wright_fisher.evolve(rng0,pop_anc,
                                  params_const_decl,recorder_const_decl) # constant pop

        fp11.wright_fisher.evolve(rng1,pop_lin_decl,
                                  params_lin_decl,recorder_lin_decl) # linear decline

        fp11.wright_fisher.evolve(rng2,pop_exp_decl,
                      params_exp_decl,recorder_exp_decl) # exponential decline

        fp11.wright_fisher.evolve(rng3,pop_disc_decl,
                      params_disc_decl,recorder_disc_decl)# discrete decline
        
        ##########
        # Use ancestral pop as the constant + outcrossing pop
        pop_const_domest_clonal = deepcopy(pop_anc) # constant + clonal prop
        pop_lin_domest_clonal = deepcopy(pop_lin_decl) # linear + outcrossing
        pop_exp_domest_clonal = deepcopy(pop_exp_decl) # exponential + outcrossing
        pop_disc_domest_clonal = deepcopy(pop_disc_decl) # discrete + outcrossing
    
    
        recorder_const_domest_out = RecordStats.RecordStats(1,rng0)    
        recorder_lin_domest_out = RecordStats.RecordStats(1,rng1) 
        recorder_exp_domest_out = RecordStats.RecordStats(1,rng2)
        recorder_disc_domest_out = RecordStats.RecordStats(1,rng3)
        recorder_const_domest_clonal = RecordStats.RecordStats(1,rng4)    
        recorder_lin_domest_clonal = RecordStats.RecordStats(1,rng5) #
        recorder_exp_domest_clonal = RecordStats.RecordStats(1,rng6)
        recorder_disc_domest_clonal = RecordStats.RecordStats(1,rng7)
        
        print("Domestication")
        fp11.wright_fisher.evolve(rng0,pop_anc,
                                  params_const_domest,recorder_const_domest_out) # constant pop  + outcrossing

        fp11.wright_fisher.evolve(rng1,pop_lin_decl,
                                  params_lin_domest,recorder_lin_domest_out) # linear decline  + outcrossing
        
        fp11.wright_fisher.evolve(rng2,pop_exp_decl,
                                  params_exp_domest,recorder_exp_domest_out) # exponential decline  + outcrossing
        
        fp11.wright_fisher.evolve(rng3,pop_disc_decl,
                                  params_disc_domest,recorder_disc_domest_out) # discrete decline  + outcrossing

        clonal.evolve(rng4,pop_const_domest_clonal,
                      params_const_domest,recorder_const_domest_clonal) # constant pop  + clonal
        
        clonal.evolve(rng5,pop_lin_domest_clonal,
                      params_lin_domest,recorder_lin_domest_clonal) # linear decline  + clonal
        
        clonal.evolve(rng6,pop_exp_domest_clonal,
                      params_exp_domest,recorder_exp_domest_clonal) # exponential decline   + clonal
        
        clonal.evolve(rng7,pop_disc_domest_clonal,
                      params_disc_domest,recorder_disc_domest_clonal) # discrete decline + clonal

        print("Writing")
        recorder_dict = {"ancestral" : recorder_anc,
                        "const_decl" : recorder_const_decl,
                        "lin_decl" : recorder_lin_decl,
                        "exp_decl" : recorder_exp_decl,
                        "disc_decl" : recorder_disc_decl,
                        "const_domest_out" : recorder_const_domest_out,
                        "lin_domest_out" : recorder_lin_domest_out,
                        "exp_domest_out" : recorder_exp_domest_out,
                        "disc_domest_out" : recorder_disc_domest_out,
                        "const_domest_clonal" : recorder_const_domest_clonal,
                        "lin_domest_clonal" : recorder_lin_domest_clonal,
                        "exp_domest_clonal" : recorder_exp_domest_clonal,
                        "disc_domest_clonal" : recorder_disc_domest_clonal}
        if (gen_model == "additive"):
            output = open('data/revisions/additive/additive_grape.'+str(rep)+'.pkl', 'wb')
        if (gen_model == "recessive"):
            output = open('data/revisions/recessive/recessive_grape.'+str(rep)+'.pkl', 'wb')
        pickle.dump(recorder_dict,output)
        output.close()
if __name__ == "__main__":
    #main()
    repid = [i for i in range(500)]
    P=mp.Pool()
    res=P.imap_unordered(main,repid)
    P.close()
    P.join()
