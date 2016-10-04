###
### Useful functions to work with PyEMMA (http://emma-project.org/latest/)
### 
### Written by Giovanni Pinamonti, 
### @SISSA, Trieste, 2015-2016


import numpy as np
import mdtraj
##########################################################################
def my_transform(self, traj):
    rad = mdtraj.compute_dihedrals(traj, self.dih_indexes, periodic=False)
    if self.cossin:
        rad = np.dstack((np.cos(rad), np.sin(rad)))
        rad = rad.reshape(rad.shape[0], rad.shape[1]*rad.shape[2])
    # convert to degrees
    if self.deg:
        rad = np.rad2deg(rad)
        
    return rad
        
#from pyemma.coordinates.data.featurizer import DihedralFeature
#DihedralFeature.transform = my_transform
##########################################################################



def get_tica_cumvar(tik):
    lastpositive=np.where(tik.eigenvalues<0)[0][0]
    totvar=np.sum(tik.eigenvalues[:lastpositive]**2)
    cumvar=[]
    for i in range(lastpositive+1):
        cumvar.append(np.sum(tik.eigenvalues[:i]**2)/totvar)
    cumvar=np.array(cumvar)
    first95=np.where(cumvar>0.95)[0][0]
    #first95=np.where(cumvar>0.9)[0][0]
    return first95,cumvar




#
# Computes the average - of a certain feature
# defined in inp - on each cluster of an MSM
#
def average_by_state(M,inp,stride=1,traj=None):
    if traj==None:
        traj=inp.get_output(stride=stride)
    #clstates=cl.index_clusters
    nfeat=inp.dimension()
    avg=np.zeros((M.nstates,nfeat))
    nn=np.zeros((M.nstates,nfeat))
    
    itraj=0
    for dt in M.discrete_trajectories_active:
        #tmp_avg=np.zeros(nfeat)
        iframe=0
        for ndx in dt[::stride]:
            avg[ndx]+=traj[itraj][iframe]
            nn[ndx]+=1
            iframe+=1
        itraj+=1
    nz=np.where(nn>0)[0] #don't divide if nn is zero
    avg[nz]/=nn[nz]      #(it can happen because of the stride)
    return avg

### This should be merged with average_by_state
### it make no sense to cmpute just the std and 
### it's slow to read the traj twice
def std_by_state(M,inp,stride=1,traj=None):
    if traj==None:
        traj=inp.get_output(stride=stride)
    #clstates=cl.index_clusters
    nfeat=inp.dimension()
    avg=np.zeros((M.nstates,nfeat))
    avg2=np.zeros((M.nstates,nfeat))
    nn=np.zeros((M.nstates,nfeat))
    
    itraj=0
    for dt in M.discrete_trajectories_active:
        #tmp_avg=np.zeros(nfeat)
        iframe=0
        for ndx in dt[::stride]:
            avg[ndx]+=traj[itraj][iframe]
            avg2[ndx]+=traj[itraj][iframe]**2
            nn[ndx]+=1
            iframe+=1
        itraj+=1
    nz=np.where(nn>0)[0] #don't divide if nn is zero
    avg[nz]/=nn[nz]      #(it can happen because of the stride)
    avg2[nz]/=nn[nz]      #(it can happen because of the stride)

    return np.sqrt(avg2-avg**2)/np.sqrt(nn) ### check here for nn zero


### In this cell I save virtual trajectories made with frames with increasing values of a TICA component
### This is useful to visualize the TIC
def print_tics(icomp,IC,inp,nframe=1000):
    print "#",icomp
    #concatenate all the TICA-projected trajectories into one
    ICtot=IC[0]
    for i in range(1,len(IC)):
        ICtot=np.concatenate((ICtot,IC[i]))
            
        #define an array with indexes ordered with growing TICA compontent icomp
    ndx=[]
    for i in range(0,len(IC)):
        for j in range(0,len(IC[i][:,icomp])):
            ndx.append(np.array([i,j]))
    ndx=np.array(ndx)
    order = np.argsort(ICtot[:,icomp])
    diff=ICtot[order[0],icomp]-ICtot[order[-1],icomp]
    #define a list of indeces equally spaced in TICA space
    minstep=abs(diff/float(nframe))
    lastframe=-999
    indeces=[]
    for i in range(len(order)):
        if lastframe+minstep>ICtot[order[i],icomp]:
            continue
        indeces.append(ndx[order[i]])
        lastframe=ICtot[order[i],icomp]
    indeces=np.array(indeces)
    print indeces.shape
    cacca=[]
    for i in range(len(indeces)):
        cacca.append(IC[indeces[i,0]][indeces[i,1],icomp])
    plt.plot(cacca,marker='.')
    pyemma.coordinates.save_traj(inp_A,indeces,'IC'+str(icomp)+'.xtc')


def dump_evec(name,M,inpXout,icomp=1,dump=False,eps=0.1,nframes=50):
### this should work fine if the fraction of active states is > 1.0 (but take a look at that)
    import matplotlib.pylab as plt
    import pyemma.coordinates as coor
    ### check with dump=False and tune the eps
    print '### HOW TO: check with dump=False and tune the eps'
    rv = M.eigenvectors_right()[:,icomp] # only active set
    n_cl=M.nstates #N. active states (good!)

    order = np.argsort(rv)
    #fig = plt.figure(figsize=(15,10))
    plt.plot(rv[order[0:n_cl]],marker='.')
    plt.plot([0,n_cl],[rv[order[0]]+eps,rv[order[0]]+eps])
    plt.plot([0,n_cl],[rv[order[-1]]-eps,rv[order[-1]]-eps])

    ndx_start=np.where(rv<rv[order[0]]+eps)[0]
    ndx_end=np.where(rv>rv[order[-1]]-eps)[0]

    fig = plt.figure(figsize=(14,3))
    plt.subplot2grid((1,2),(0,0))
    plt.plot(rv[ndx_start], marker='o',lw=0.1)
    plt.title('START')
    plt.subplot2grid((1,2),(0,1))
    plt.plot(rv[ndx_end], marker='o',lw=0.1)
    plt.title('END')
    
    start=np.zeros(n_cl)
    start[ndx_start]+=1
    start=start/len(ndx_start)
    end=np.zeros(n_cl)
    end[ndx_end]+=1
    end=end/len(ndx_end)
#    print ndx_start
#    print start
#    return
    start_frames=M.sample_by_distributions([start],nframes)[0]
    end_frames=M.sample_by_distributions([end],nframes)[0]
#    start=cl.sample_indexes_by_cluster(ndx_start,nframes)
#    end=cl.sample_indexes_by_cluster(ndx_end,nframes)
#    print start_frames
    print 'Number of microstates: start ->',len(start_frames),'  end ->', len(end_frames)
    if dump:
        coor.save_traj(inpXout,start_frames,name+'_evec'+str(icomp)+'-start.pdb')
        coor.save_traj(inpXout,end_frames,name+'_evec'+str(icomp)+'-end.pdb')
        
    return


def dump_process(name,M,inpXout,icomp=1,nframes=10):
### this should work fine if the fraction of active states is > 1.0 (but take a look at that)
    import matplotlib.pylab as plt
    import pyemma.coordinates as coor
    ### check with dump=False and tune the eps
    print '### HOW TO: check with dump=False and tune the eps'
    rv = M.eigenvectors_right()[:,icomp] # only active set
    n_cl=M.nstates #N. active states (good!)

    order = np.argsort(rv)
    #fig = plt.figure(figsize=(15,10))
    plt.plot(rv[order[0:n_cl]],marker='.')
    
    frames=np.array([])
    for icl in order:
        frames=np.concatenate(frames,M.sample_by_distributions([icl],nframes)[0])

    print 'Number of frames dumped =', len(frames)

    coor.save_traj(inpXout,frames,name+'_proc'+str(icomp)+'.pdb')
        
    return


def plot_3traj(M,traj_x,traj_y,traj_z,nbins=50):
    #### Plot using traj_x, traj_y as x,y coordinates,
    #### traj_z as color (averaged over each bin)
    F=np.zeros((nbins,nbins))
    NN=np.zeros((nbins,nbins))
    
    
    extent=np.zeros(4)
    extent[0]= np.min(traj_x)#min x
    extent[1]= np.max(traj_x)#max x
    extent[2]= np.min(traj_y)#min y
    extent[3]= np.max(traj_y)#max y
    dx=(extent[1]-extent[0])/nbins
    dy=(extent[3]-extent[2])/nbins
    k=0
    for g,d in zip(traj_x,traj_y):
        ix=int((g-extent[0])/dx)
        if ix>=nbins: ix=nbins-1
        iy=int((d-extent[2])/dy)
        if iy>=nbins: iy=nbins-1
        F[ix,iy]+=traj_z[k]
        NN[ix,iy]+=1
        k+=1
    F/=NN
    return F,extent



def _evaluate_predictions(T,hmm_obj,nsteps,stride=1):
### In this cell I evaluate the probality to go from set A to set B 
### at time k*tau as predicted by a certain MSM transition matrix
    nstates=hmm_obj.nstates
    size=T.shape[0]
    TK=np.identity(size)
    PK=[]
    for k in range(nsteps):   
        TK=np.matmul(TK,T)
        myP=np.zeros((nstates,nstates))
        if not k%stride==0: continue
        for A in range(nstates):
            for B in range(nstates):
                for j in range(size):
                    tmp_P=np.dot(TK[:,j],hmm_obj.metastable_distributions[A,:])
                    tmp_P*=hmm_obj.metastable_memberships[j,B]
                    myP[A,B]+=tmp_P
        PK.append(myP)
    PK=np.array(PK)
    return PK

def predicted_transitions(M,hmm_obj,nsteps,nsamples=50,stride=1):
### PREDICTED TRANSITIONS
### Here I sample 10 transition matrices
### then I propagate them for 100 steps and
### I store the probabilities of transition
### between metastable sets (of M)
    ctrajs=M.discrete_trajectories_full
    lag=M.lagtime
    import msmtools
    C=msmtools.estimation.effective_count_matrix(ctrajs,1)
    T0=msmtools.estimation.tmatrix(msmtools.estimation.connected_cmatrix(C),reversible=True)
    sampleMs=msmtools.estimation.tmatrix_sampler(\
                msmtools.estimation.connected_cmatrix(C),reversible=True,T0=T0.toarray())
    sampledPKest=[]
    for i in range(nsamples):
        T=sampleMs.sample()
        tmp_PK=_evaluate_predictions(T,hmm_obj,nsteps,stride=stride)
        sampledPKest.append(tmp_PK)
    return np.array(sampledPKest)


def estimated_transitions(M,hmm_obj,nsteps,nsamples=50,stride=1):
### ESTIMATED TRANSITIONS
### Here, for lagtimes multiples of  
### M.lagtime, I sample 10 MSM and
### I store the probabilities of transition
### between metastable sets (of M)
    ctrajs=M.discrete_trajectories_full
    lag=M.lagtime
    import msmtools
    sampledEKest=[]
    nstates=hmm_obj.nstates
    for i in range(1,nsteps+1,stride):
        CK=msmtools.estimation.effective_count_matrix(ctrajs,lag*i) #this is slow (56 s)
        TK0=msmtools.estimation.transition_matrix(\
                    msmtools.estimation.connected_cmatrix(CK),reversible=True)
        sampleMKs=msmtools.estimation.tmatrix_sampler(
            msmtools.estimation.connected_cmatrix(CK),reversible=True,T0=TK0.toarray())
        tmp_EK=[]
        active_set=msmtools.estimation.connected_sets(CK)[0]
        size=active_set.shape[0]
        for i in range(nsamples):
            TK=sampleMKs.sample()
            myP=np.zeros((nstates,nstates))
            for A in range(nstates):
                for B in range(nstates):
                    for j in range(size):
                        tmp_P=np.dot(TK[:,j],\
                                hmm_obj.metastable_distributions[A,active_set]\
                                /np.sum(hmm_obj.metastable_distributions[A,active_set]))
                        tmp_P*=hmm_obj.metastable_memberships[active_set[j],B]
                        myP[A,B]+=tmp_P
            tmp_EK.append(myP)
        sampledEKest.append(tmp_EK)
    return np.array(sampledEKest)
