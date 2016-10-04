# ------------------------------------------------------------
# This functions are taken from the software baRNAba (https://github.com/srnas/barnaba)
# They compute the G-vectors and are compatible with
# pyemma.coordinates.featurizer (see http://emma-project.org/latest/ for PyEMMA documentation)
#
# Written by Giovanni Pinamonti, Berlin, 29 Ottobre 2015
#
# G-vectors have been defined in
# S.Bottaro, F. di Palma, G.Bussi. The role of nucleobases in RNA structure and dynamics. Nucleic Acids Research (2014)
# Please cite it if you use this code
#
# NB: baRNAba reads PDB and uses angstrom as unit of measure
#     while mdtraj uses nanometers. So I rescaled the cutoff radius
#     from 2.4 to 0.24
# ------------------------------------------------------------

import mdtraj
from scipy.spatial import distance
import numpy as np
import sys

f_factors=[5.,5.,3.]
scale=[1./f_factors[0],1./f_factors[1],1./f_factors[2]]
R_C=0.24
#R_C=0.32

def get_lcs_idx(top):
    pur_list=['A','G']
    pyr_list=['U','C']
#    print top.n_residues
    CA_idx = []
    CB_idx = []
    CC_idx = []
    for res in top.residues:
        i1=[idx.index for idx in res.atoms_by_name('C2')]
        #assert len(i1) == 1, "# missing C2 in %s \n" % (self.sequence_id[i])
        assert len(i1) == 1, "# missing C2 \n"
        i1=i1[0]
        
        i2=[idx.index for idx in res.atoms_by_name('C4')]
        #assert len(i2) == 1, "# missing C4 in %s \n" % (self.sequence_id[i])
        assert len(i2) == 1, "# missing C4 \n"
        i2=i2[0]
        
        i3=[idx.index for idx in res.atoms_by_name('C6')]
        #assert len(i3) == 1, "# missing C6 in %s \n" % (self.sequence_id[i])
        assert len(i3) == 1, "# missing C6 \n"
        i3=i3[0]
        
        CA_idx.append(i1)
        
        mytype = res.name
        if mytype in pyr_list:
            CB_idx.append(i2)
            CC_idx.append(i3)
        if mytype in pur_list:
            CB_idx.append(i3)
            CC_idx.append(i2)

    lcs_idx =  np.array([CA_idx,CB_idx, CC_idx])
    return lcs_idx

def get_lcs(lcs_idx,coords):

    coords1 = coords[lcs_idx[0]]
    coords2 = coords[lcs_idx[1]]
    coords3 = coords[lcs_idx[2]]
    coords = np.array([coords1,coords2,coords3])

    # calculate center of mass
    origo = np.sum(coords,axis=0)/3.0
        
    # CoM-C2 (x axis)
    x = coords[0]-origo
    x_norm = np.sqrt(np.sum(x*x,axis=1))
    x = x/x_norm[:,np.newaxis]
    # CoM-C4/C6 
    c = coords[1]-origo
    # z/y axis
    z = np.cross(x,c,axis=1)
    z_norm = np.sqrt(np.sum(z*z,axis=1))
    z = z/z_norm[:,np.newaxis]
    y = np.cross(z,x,axis=1)
    lcs = np.array([x.T,y.T,z.T]).T
    return lcs,origo

def get_3dmat(lcs_idx,coords,cutoff):

    lcs, origo = get_lcs(lcs_idx,coords)

    # prune search first
    max_r  = np.max(f_factors)*cutoff ### questo si puo' fare fuori
    dmat = distance.squareform(distance.pdist(origo))
    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.001))).T
#    m_idx = np.array(np.where((dmat<max_r) & (dmat>0.01))).T
        
    # calculate scaled distances
    diff = [origo[y]-origo[x] for x,y in m_idx]
    dotp = np.array([np.dot(diff[i],lcs[j]) for i,j in zip(range(len(diff)),m_idx[:,0])])

    return dotp,m_idx

def get_gmat(lcs_idx,coords,cutoff):
#    print lcs_idx.shape
    ll=lcs_idx.shape[1]
    mat = np.zeros((ll,ll,4))
            
    dotp,m_idx = get_3dmat(lcs_idx,coords,cutoff)
        # if all zeros
    if(len(dotp)==0):
        return mat

    dotp *= np.array(scale)[np.newaxis,:]

    dotp_norm = np.sqrt(np.sum(dotp**2,axis=1))
        
    # calculate 4D g-vector
    ff = (np.pi*dotp_norm)/cutoff
    factor13 = np.sin(ff)/ff
    factor4= ((1.0+np.cos(ff))*cutoff)/np.pi
    gmat = dotp*factor13[:,np.newaxis]
    gmat = np.concatenate((gmat,factor4[:,np.newaxis]),axis=1)

    # set to zero when norm is larger than cutoff
    gmat[dotp_norm>cutoff] = 0.0
        
    #mat = np.zeros((ll,ll,4))
    mat[m_idx[:,0],m_idx[:,1]] = gmat
        
    return mat


# function get_gvecs
# to be used in pyemma.coordinates.featurizer.add_custom_func(get_gvecs,4*nres*nres)
# input: mdtraj.trajectory object
# output: numpy.array (compatible with pyemma featurizer)
def get_gvecs(traj):
    top=traj.topology
    l_nb=top.select('name C2 or name C4 or name C6')
    t_nb=traj.atom_slice(l_nb)

    nres=t_nb.top.n_residues
    assert nres > 1, "# %d residue selected. Cannot compute G-vectors \n" % (nres)

    lcs_idx=get_lcs_idx(t_nb.top)
    gvecs=[]

    for iframe in range(len(traj)):
        coords=t_nb.slice(iframe).xyz[0]
#        gmat=get_gmat(lcs_idx,coords,0.24)
        gmat=get_gmat(lcs_idx,coords,R_C)
        gvecs.append(gmat)
        
    gvecs=np.array(gvecs,dtype=np.float32).reshape(len(traj),4*nres*nres)
    return gvecs

# function get_gvecs with non standard cutoff radius
def get_gvecs_r32(traj,cutoff):
    top=traj.topology
    l_nb=top.select('name C2 or name C4 or name C6')
    t_nb=traj.atom_slice(l_nb)

    nres=t_nb.top.n_residues
    assert nres > 1, "# %d residue selected. Cannot compute G-vectors \n" % (nres)

    lcs_idx=get_lcs_idx(t_nb.top)
    gvecs=[]

    for iframe in range(len(traj)):
        coords=t_nb.slice(iframe).xyz[0]
        gmat=get_gmat(lcs_idx,coords,cutoff) ### DIFFERENT CUTOFF RADIUS
        gvecs.append(gmat)
        
    gvecs=np.array(gvecs,dtype=np.float32).reshape(len(traj),4*nres*nres)
    return gvecs



def get_sel_gvecs(traj,stringa):
    if len(stringa)==0:
        return get_gvecs(traj)
    else:
        l_res=top.select(stringa)
        tr_res=traj.atom_slice(l_res)
        tp_res=tr_res.topology
        l_nb=tp_res.select('name C2 or name C4 or name C6')
        t_nb=tr_res.atom_slice(l_nb)
        return get_gvecs(t_nb)    
