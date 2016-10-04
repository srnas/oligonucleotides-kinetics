# ------------------------------------------------------------
# Functions to quickly read/write pyemma objects
#
# Written by Giovanni Pinamonti 
#  SISSA, Trieste, Italy, 2016
#
# ------------------------------------------------------------



import numpy as np
import cPickle as pickle

###TICA write
def write_tica_eval(name,tica_obj):
    fh=open(name+'_tica_eval.dat','w')
    fh.write("# TICA with lagtime="+str(tica_obj.lag)+" steps\n")
    for eval in tica_obj.eigenvalues:
        fh.write("%e\n" %(eval))
    fh.close()

def write_tica_evec(name,tica_obj):
    fh=open(name+'_tica_evec.dat','w')
    fh.write("# TICA with lagtime="+str(tica_obj.lag)+" steps\n")
    for evec in tica_obj.eigenvectors:
        for x in evec:
            fh.write("%e " %(x))
        fh.write("\n")
    fh.close()

def write_tica_trajs(name,IC):
    itrj=0
    for tictrj in IC:
        fh=open(name+'_tica_traj'+str(itrj)+'.dat','w')
#        fh.write("# TICA with lagtime="+str(tica_obj.lag)+" steps\n")
        for frame in tictrj:
            for x in frame:
                fh.write(str(x)+" ")
            fh.write("\n")
        fh.close()
        itrj+=1

###TICA read
def read_tica_eval(name):
    ### to read eigenvalues use this
    tica_eval=[]
    for line in open(name+'_tica_eval.dat','r'):
        if line[0]=='#':
            continue
        tica_eval.append(np.array(float(line.split()[0])))
    return np.array(tica_eval)

def read_tica_evec(name):
    ### to read eigenvectors use this
    tica_evec=[]
    for line in open(name+'_tica_evec.dat','r'):
        if line[0]=='#':
            continue
        tica_evec.append(np.array([float(x) for x in line.split()]))
    return np.array(tica_evec)

def read_tica_trajs(name,ntraj):
    ### to read TIC-trajectories use this
    NTRAJ=ntraj
    IC=[]
    for itrj in range(NTRAJ):
        tica_trj=[]
        for line in open(name+'_tica_traj'+str(itrj)+'.dat','r'):
            if line[0]=='#':
                continue
            tica_trj.append(np.array([float(x) for x in line.split()]))
        tica_trj=np.array(tica_trj)
        IC.append(tica_trj)
    return IC



#######################################
############ CLUSTERING ###############
#######################################

def write_cl_dtrajs(name,dtrajs):
    itrj=0
    for trj in dtrajs:
        fh=open(name+'_cl_dtraj'+str(itrj)+'.dat','w')
        for x in trj:
            fh.write(str(x)+"\n")
        fh.close()
        itrj+=1

def write_cl_centers(name,cl):
    fh=open(name+'_cl_centers.dat','w')
    fh.write("# Clustering details: "+cl.describe()+"\n")
    for cc in cl.clustercenters:
        for x in cc:
            fh.write(str(x)+" ")
        fh.write("\n")
    fh.close()

import cPickle as pickle
def write_cl_indexes(name,cl):
    fh=open(name+'_cl_indexes.idx','w')
    pickle.dump(cl.index_clusters,fh,-1)
    fh.close()

def read_cl_indexes(name):
    fh=open(name+'_cl_indexes.idx','r')
    cl_indexes=pickle.load(fh)
    fh.close()
    return cl_indexes


def read_cl_dtrajs(name,ntraj):
    ### to read discrete trajectories use this
    NTRAJ=ntraj
    dtrajs=[]
    for itrj in range(NTRAJ):
        trj=[]
        for line in open(name+'_cl_dtraj'+str(itrj)+'.dat','r'):
            if line[0]=='#':
                continue
            trj.append(int(line))
        trj=np.array(trj)
        dtrajs.append(trj)
    return dtrajs

def read_cl_centers(name):
    ### to read cluster centers use this
    cl_centers=[]
    for line in open(name+'_cl_centers.dat','r'):
        if line[0]=='#':
            continue
        cl_centers.append(np.array([float(x) for x in line.split()]))
    return np.array(cl_centers)




######################################
############   M.S.M.  ###############
######################################
def write_msm_tmat(name,M):
    fh=open(name+'_msm_tmat.dat','w')
    fh.write("# MSM with a lagtime of "+str(M.lag)+" steps\n")
    for line in M.transition_matrix:
        for x in line:
            fh.write(str(x)+" ")
        fh.write("\n")
    fh.close()

def read_msm_tmat(name):
    ### to read msm's transition matrix
    T=[]
    for line in open(name+'_msm_tmat.dat','r'):
        if line[0]=='#':
            continue
        T.append(np.array([float(x) for x in line.split()]))
    T=np.array(T)
    if T.shape[0]!=T.shape[1]:
        print "WARNING: T is not q square matrix!"
    return T



##########################################333
##############3 timescales!!!!!!1 #########33
#####################3333#################333

def write_msm_its(name,its):
    fh=open(name+'_msm_its.dat','w')
    for lag,ts in zip(its.lags,its.timescales):
        fh.write(str(lag)+' ')
        for x in ts:
            fh.write(str(x)+' ')
        fh.write('\n')
    fh.close()

def write_hmm_its(name,hmits):
    fh=open(name+'_hmm_its.dat','w')
    for lag,ts in zip(hmits.lags,hmits.timescales):
    #print lag,ts
        fh.write(str(lag)+' ')
        for x in ts:
            fh.write(str(x)+' ')
        fh.write('\n')
    fh.close()

def read_msm_its(name):
    lagtimes=[]
    timescales=[]
    for line in open(name+'_msm_its.dat','r'):
        if line[0]=='#':
            continue
        lag=int(line.split()[0])
        ts=[float(x) for x in line.split()[1:]]
        lagtimes.append(lag)
        timescales.append(np.array(ts))
        
    return np.array(lagtimes),np.array(timescales)


def read_hmm_its(name):
    lagtimes=[]
    timescales=[]
    for line in open(name+'_hmm_its.dat','r'):
        if line[0]=='#':
            continue
        lag=int(line.split()[0])
        ts=[float(x) for x in line.split()[1:]]
        lagtimes.append(lag)
        timescales.append(np.array(ts))
        
    return np.array(lagtimes),np.array(timescales)

def write_its_errors(name,its):
    fh=open(name+'_msm_its_errors.dat','w')
    fh.write('# lag, avg, std')
    for lag,ts,er in zip(its.lags,its.sample_mean,its.sample_std):
        fh.write(str(lag)+' ')
        for x in ts:
            fh.write(str(x)+' ')
        for x in er:
            fh.write(str(x)+' ')
        fh.write('\n')
    fh.close()

def write_hmm(name,HiddenMM):
    fh=open(name+'.hmm','wb')
    pickle.dump(HiddenMM,fh,-1)
    fh.close()   

def read_hmm(name):
    return pickle.load(open(name+'.hmm','rb'))

def write_msm(name,M):
    fh=open(name+'.msm','wb')
    pickle.dump(M,fh,-1)
    fh.close()   

def read_msm(name):
    return pickle.load(open(name+'.msm','rb'))


def write_dpa_cl(name,cl_dpa):
    fh=open(name+'.dpc','wb')
    pickle.dump(cl_dpa,fh,-1)
    fh.close()   

def read_dpa_cl(name):
    return pickle.load(open(name+'.dpc','rb'))
