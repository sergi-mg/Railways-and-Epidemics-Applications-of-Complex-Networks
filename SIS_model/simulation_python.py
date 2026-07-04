#code for simulating
#%%
#import libraries
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from numba import njit
import networkx as nx
import random as random
import subprocess
import sys
import time
import os
from scipy.optimize import curve_fit


#%%
#functions
def sample_powerlaw(N, gamma, kmin, kc):
    k = np.arange(kmin, kc + 1, 1, dtype=int)
    p = k ** (-gamma)
    p /= p.sum()

    return np.random.choice(k, size=N, p=p)

def G_to_PD(G):

    #change of saving structure: vector and pointer
    V=np.array(list(dict(G.degree()).values()),dtype=int)
    P=np.zeros((V.size+1),dtype=int)

    P[1:]=np.cumsum(V)
    
    D = np.empty(P[-1], dtype=int)
    pair = np.empty(P[-1], dtype=int)


    edge_id = {}  # (u,v) -> index

    idx = 0

    for u in range(N):
        neighbours=G.adj[u].keys()
        for v in neighbours:

            D[idx] = v
            edge_id[(u, v)] = idx
            idx += 1

    # construir inversos
    for (u, v), i in edge_id.items():
        pair[i] = edge_id[(v, u)]
        

    return P,D,pair


# %%
#parameters
N_list=np.array([10000, 30000, 50000,100000, 300000, 500000, 1000000])
N_list=N_list[0:0]
N_sim=2000
N_networks=5

#%%
gamma=2.5
lambda_min=np.array([0.032,0.024,0.016,0.012,0.012,0.008,0.007])
lambda_max=np.array([0.035,0.026,0.020,0.017,0.013,0.010,0.009])
lambda_min=lambda_min[0:0]
lambda_max=lambda_max[0:0]
d_lambda=0.0002

#simulations

for N_i in range(len(N_list)):
    N=N_list[N_i]
    k_c=(N)**(1/(gamma-1))
    lambda_list=np.arange(lambda_min[N_i],lambda_max[N_i]+0.0001,d_lambda)
    N_lambda=lambda_list.size
    data_results=np.zeros((N_sim*N_networks,N_lambda,2),dtype=float)
    
    for i in range(N_networks):
        position=i*N_sim
        #create the network
        np.random.seed(i)
        random.seed(i)

        #vector of k values following a power law
        vector_k=sample_powerlaw(N,gamma,4,k_c)
        E_2=np.sum(vector_k)
        if E_2%2!=0:
            vector_k[np.random.randint(N)]+=1

        #configuration model
        G = nx.configuration_model(vector_k,seed=i)
        G = nx.Graph(G)
        G.remove_edges_from(nx.selfloop_edges(G))

        #change in format
        P,D,D_sim=G_to_PD(G)
        P+=1 
        D+=1
        D_sim+=1
        E=D.size//2

        #save to files
        fD_sim = open("../data_networks/D_sim"+str(i)+".dat", "w")
        np.savetxt(fD_sim, D_sim, fmt="%d")           
        fD_sim.flush()
        os.fsync(fD_sim.fileno())
        fD_sim.close()
        
        fP = open("../data_networks/P"+str(i)+".dat", "w")
        np.savetxt(fP, P, fmt="%d")           
        fP.flush()
        os.fsync(fP.fileno())
        fP.close()
        
        fD = open("../data_networks/D"+str(i)+".dat", "w")
        np.savetxt(fD, D, fmt="%d")
        fD.flush()
        os.fsync(fD.fileno())
        fD.close()  
        
        for j in range(N_lambda):
            lam=lambda_list[j]
        
            cmd = ["./net2.exe", str(N), str(gamma)+"d0",
                   str(lam)+"d0", str(i), str(E), str(N_sim)]
            
            result = subprocess.run(cmd)
            
            results=np.loadtxt("../data_results_fortran/results_simulation.dat"
                               ,usecols=[-2,-1])
            
            data_results[position:position+N_sim,j,:]=results
    
    for j in range(N_lambda): 
        lam=lambda_list[j]

        #save the results for each (N,lambda)
        name="../data_results/N_"+str(N)+"_gamma_"+str(round(gamma,1))\
        +"_lambda_"+str(round(lam,4))+"_Nsim_"+str(N_sim)+".dat" 
            
        np.savetxt(name, data_results[:,j,:])  
            
   
#%%

N_list=np.array([10000, 30000, 50000,100000, 300000, 500000, 1000000])
N_list=N_list[[3,4]]
gamma=3.5
lambda_min=np.array([0.15,0.132,0.120,0.122,0.114,0.106,0.101])
lambda_max=np.array([0.16,0.14,0.130,0.130,0.115,0.110,0.105])
lambda_min=lambda_min[[3,4]]
lambda_max=lambda_max[[3,4]]
d_lambda=0.0002

#simulations

for N_i in range(len(N_list)):
    N=N_list[N_i]
    k_c=(N)**(1/(gamma-1))
    lambda_list=np.arange(lambda_min[N_i],lambda_max[N_i]+0.0001,d_lambda)
    N_lambda=lambda_list.size
    
    data_results=np.zeros((N_sim*N_networks,N_lambda,2),dtype=float)
    
    for i in range(N_networks):
        position=i*N_sim
        #create the network
        np.random.seed(i)
        random.seed(i)

        #vector of k values following a power law
        vector_k=sample_powerlaw(N,gamma,4,k_c)
        E_2=np.sum(vector_k)
        if E_2%2!=0:
            vector_k[np.random.randint(N)]+=1

        #configuration model
        G = nx.configuration_model(vector_k,seed=i)
        G = nx.Graph(G)
        G.remove_edges_from(nx.selfloop_edges(G))
        
        #change in format
        P,D,D_sim=G_to_PD(G)
        P+=1 
        D+=1
        D_sim+=1
        E=D.size//2

        #save to files
        fD_sim = open("../data_networks/D_sim"+str(i)+".dat", "w")
        np.savetxt(fD_sim, D_sim, fmt="%d")           
        fD_sim.flush()
        os.fsync(fD_sim.fileno())
        fD_sim.close()
        
        fP = open("../data_networks/P"+str(i)+".dat", "w")
        np.savetxt(fP, P, fmt="%d")           
        fP.flush()
        os.fsync(fP.fileno())
        fP.close()
        
        fD = open("../data_networks/D"+str(i)+".dat", "w")
        np.savetxt(fD, D, fmt="%d")
        fD.flush()
        os.fsync(fD.fileno())
        fD.close()  
        
        for j in range(N_lambda):
            lam=lambda_list[j]
        
            cmd = ["./net2.exe", str(N), str(gamma)+"d0",
                   str(lam)+"d0", str(i), str(E), str(N_sim)]
            
            result = subprocess.run(cmd)
            
            results=np.loadtxt("../data_results_fortran/results_simulation.dat"
                               ,usecols=[-2,-1])
            
            data_results[position:position+N_sim,j,:]=results
    
    for j in range(N_lambda): 
        lam=lambda_list[j]

        #save the results for each (N,lambda)
        name="../data_results/N_"+str(N)+"_gamma_"+str(round(gamma,1))\
        +"_lambda_"+str(round(lam,4))+"_Nsim_"+str(N_sim)+".dat" 
            
        np.savetxt(name, data_results[:,j,:])  
        
