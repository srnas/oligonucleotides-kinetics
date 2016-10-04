### ------------------------------
### Definitions of stacking score 
### (see Condon et al., JCTC 2014, 11, 2729-2742)
### compatibles with pyemma featurizer
### written by Giovanni Pinamonti (SISSA, Trieste)
### @FU-Berlin, November 2015
### ------------------------------

import numpy as np
import sys
import mdtraj

atom_name_base='name N1 or name C2 or name N2 or name O2 or name N3 or name C4 or name N4 or name O4 or name C5 or name C6 or name N6 or name O6 or name N7 or name C8 or name N9'

def get_comdist(traj,res1,res2):
    top=traj.topology
    l_b1=top.select('resid '+str(res1)+' and ('+atom_name_base+')')
    l_b2=top.select('resid '+str(res2)+' and ('+atom_name_base+')')

    tb1=traj.atom_slice(l_b1)
    tb2=traj.atom_slice(l_b2)
    pos1=tb1.xyz
    pos2=tb2.xyz
    d0=[]
    natoms1=len(l_b1)
    natoms2=len(l_b2)
    for iframe in range(len(pos1)):
        tmp_CoM1=np.array([sum(pos1[iframe,:,0]),sum(pos1[iframe,:,1]),sum(pos1[iframe,:,2])])/natoms1
        tmp_CoM2=np.array([sum(pos2[iframe,:,0]),sum(pos2[iframe,:,1]),sum(pos2[iframe,:,2])])/natoms2
        vec_d0=tmp_CoM1-tmp_CoM2
        tmp_d0=np.sqrt(sum(vec_d0**2))
        d0.append(tmp_d0)
        
    d0=np.array(d0,dtype=np.float32).reshape(len(traj),1)
    return d0

def get_omega21(traj,res1,res2):
    top=traj.topology
    l_nb1=top.select('resid '+str(res1)+' and ('+atom_name_base+')')
    l_nb2=top.select('resid '+str(res2)+' and ('+atom_name_base+')')
 
    tnb1=traj.atom_slice(l_nb1)
    tnb2=traj.atom_slice(l_nb2) 
    #resname2=tnb2.top.residue(res2).name
    resname2=top.residue(0).name
    if resname2=='A':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name N6')
    if resname2=='G':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name O6')
    if resname2=='C':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name N4')
    if resname2=='U':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name O4')

        
    pos1=tnb1.xyz
    pos2=tnb2.xyz
    omegas=[]
    natoms1=tnb1.n_atoms
    natoms2=tnb2.n_atoms
    
    a2s=tnb2.atom_slice(l_C8_2).xyz
    b2s=tnb2.atom_slice(l_N6_2).xyz
    for iframe in range(len(traj)):
        CoM1=np.array([sum(pos1[iframe,:,0]),sum(pos1[iframe,:,1]),sum(pos1[iframe,:,2])])/natoms1
        CoM2=np.array([sum(pos2[iframe,:,0]),sum(pos2[iframe,:,1]),sum(pos2[iframe,:,2])])/natoms2
        d0=CoM1-CoM2
        a2=a2s[iframe,0,:]-CoM2
        b2=b2s[iframe,0,:]-CoM2
        c2=np.cross(a2,b2)
        len_c2=np.sqrt(sum(c2**2))
        len_d0=np.sqrt(sum(d0**2))
        omega=np.arccos(np.dot(c2,d0)/(len_c2*len_d0))
        omegas.append(omega)
    omegas=np.array(omegas,dtype=np.float32).reshape(len(traj),1)
    return omegas

def get_omega12(traj,res1,res2):
    top=traj.topology
    l_nb1=top.select('resid '+str(res1)+' and ('+atom_name_base+')')
    l_nb2=top.select('resid '+str(res2)+' and ('+atom_name_base+')')

    tnb1=traj.atom_slice(l_nb1)
    tnb2=traj.atom_slice(l_nb2)    
    resname1=tnb1.top.residue(0).name
    if resname1=='A':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name N6')
    if resname1=='G':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name O6')
    if resname1=='C':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name N4')
    if resname1=='U':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name O4')
        
    pos1=tnb1.xyz
    pos2=tnb2.xyz
    omegas=[]
    natoms1=tnb1.n_atoms
    natoms2=tnb2.n_atoms
    a1s=tnb1.atom_slice(l_C8_1).xyz
    b1s=tnb1.atom_slice(l_N6_1).xyz
    for iframe in range(len(traj)):
        CoM1=np.array([sum(pos1[iframe,:,0]),sum(pos1[iframe,:,1]),sum(pos1[iframe,:,2])])/natoms1
        CoM2=np.array([sum(pos2[iframe,:,0]),sum(pos2[iframe,:,1]),sum(pos2[iframe,:,2])])/natoms2
        d0=CoM2-CoM1
        a1=a1s[iframe,0,:]-CoM1
        b1=b1s[iframe,0,:]-CoM1
        c1=np.cross(a1,b1)
        len_c1=np.sqrt(sum(c1**2))
        len_d0=np.sqrt(sum(d0**2))
        omega=np.arccos(np.dot(c1,d0)/(len_c1*len_d0))
        #omegas.append(len_c1)
        omegas.append(omega)
    omegas=np.array(omegas,dtype=np.float32).reshape(len(traj),1)
    return omegas

def get_chi(traj,res1,res2):
    top=traj.topology
    l_nb1=top.select('resid '+str(res1)+' and ('+atom_name_base+')')
    l_nb2=top.select('resid '+str(res2)+' and ('+atom_name_base+')')

    tnb1=traj.atom_slice(l_nb1)
    tnb2=traj.atom_slice(l_nb2)    
        
    resname1=tnb1.top.residue(0).name
    if resname1=='A':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name N6')
    if resname1=='G':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name O6')
    if resname1=='C':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name N4')
    if resname1=='U':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name O4')
    resname2=tnb2.top.residue(0).name
    if resname2=='A':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name N6')
    if resname2=='G':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name O6')
    if resname2=='C':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name N4')
    if resname2=='U':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name O4')

    pos1=tnb1.xyz
    pos2=tnb2.xyz
    chis=[]
    natoms1=tnb1.n_atoms
    natoms2=tnb2.n_atoms
    a1s=tnb1.atom_slice(l_C8_1).xyz
    b1s=tnb1.atom_slice(l_N6_1).xyz
    a2s=tnb2.atom_slice(l_C8_2).xyz
    b2s=tnb2.atom_slice(l_N6_2).xyz
    for iframe in range(len(traj)):
        CoM1=np.array([sum(pos1[iframe,:,0]),sum(pos1[iframe,:,1]),sum(pos1[iframe,:,2])])/natoms1
        CoM2=np.array([sum(pos2[iframe,:,0]),sum(pos2[iframe,:,1]),sum(pos2[iframe,:,2])])/natoms2
        a1=a1s[iframe,0,:]-CoM1
        b1=b1s[iframe,0,:]-CoM1
        c1=np.cross(a1,b1)
        len_c1=np.sqrt(sum(c1**2))
        a2=a2s[iframe,0,:]-CoM2
        b2=b2s[iframe,0,:]-CoM2
        c2=np.cross(a2,b2)
        len_c2=np.sqrt(sum(c2**2))
        chi=np.arccos(np.dot(c1,c2)/(len_c1*len_c2))
        chis.append(chi)
    chis=np.array(chis,dtype=np.float32).reshape(len(traj),1)
    return chis

def get_stack_score(traj,res1,res2):
    top=traj.topology
    l_nb1=top.select('resid '+str(res1)+' and ('+atom_name_base+')')
    l_nb2=top.select('resid '+str(res2)+' and ('+atom_name_base+')')

    tnb1=traj.atom_slice(l_nb1)
    tnb2=traj.atom_slice(l_nb2)
    pos1=tnb1.xyz
    pos2=tnb2.xyz
    score=[]
    natoms1=tnb1.n_atoms
    natoms2=tnb2.n_atoms
    
    resname1=tnb1.top.residue(0).name
    if resname1=='A':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name N6')
    if resname1=='G':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name O6')
    if resname1=='C':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name N4')
    if resname1=='U':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name O4')
    resname2=tnb2.top.residue(0).name
    if resname2=='A':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name N6')
    if resname2=='G':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name O6')
    if resname2=='C':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name N4')
    if resname2=='U':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name O4')

    a2s=tnb2.atom_slice(l_C8_2).xyz
    b2s=tnb2.atom_slice(l_N6_2).xyz
    a1s=tnb1.atom_slice(l_C8_1).xyz
    b1s=tnb1.atom_slice(l_N6_1).xyz
 
    for iframe in range(len(pos1)):
        #compute d0
        CoM1=np.array([sum(pos1[iframe,:,0]),sum(pos1[iframe,:,1]),sum(pos1[iframe,:,2])])/natoms1
        CoM2=np.array([sum(pos2[iframe,:,0]),sum(pos2[iframe,:,1]),sum(pos2[iframe,:,2])])/natoms2
        d0=CoM1-CoM2
        len_d0=np.sqrt(sum(d0**2))
        #compute omega21
        a2=a2s[iframe,0,:]-CoM2
        b2=b2s[iframe,0,:]-CoM2
        c2=np.cross(a2,b2)
        len_c2=np.sqrt(sum(c2**2))
        omega=np.arccos(np.dot(c2,d0)/(len_c2*len_d0))*180/np.pi
        #compute chi
        a1=a1s[iframe,0,:]-CoM1
        b1=b1s[iframe,0,:]-CoM1
        c1=np.cross(a1,b1)
        len_c1=np.sqrt(sum(c1**2))
        chi=np.arccos(np.dot(c1,c2)/(len_c1*len_c2))*180/np.pi
        #now compute the score:
        #d0
        if len_d0<=0.35:
            tmp_score=1
        elif len_d0 <=0.5:
            #tmp_score=(0.5-len_d0)/0.15
            #tmp_score=(0.5-len_d0)**3/0.15**3
            A=0.065258752
            tmp_score=A/(len_d0**3)-8*A
        else:
            score.append(0)
            continue
        #omega
        if omega<=25 or omega>=155:
            tmp_score+=1
        elif omega<=50:
            tmp_score+=(50-omega)/25
        elif omega>=130:
            tmp_score+=(omega-130)/25
        else:
            score.append(0)
            continue
        #chi
        if chi>45 and chi<135:
            tmp_score*=-1        
        #append final score
        score.append(tmp_score)
        
    score=np.array(score,dtype=np.float32).reshape(len(traj),1)
    return score


def get_stack_var(traj,res1,res2):
    top=traj.topology
    l_nb1=top.select('resid '+str(res1)+' and ('+atom_name_base+')')
    l_nb2=top.select('resid '+str(res2)+' and ('+atom_name_base+')')

    tnb1=traj.atom_slice(l_nb1)
    tnb2=traj.atom_slice(l_nb2)
    pos1=tnb1.xyz
    pos2=tnb2.xyz
    varis=[]
#    len_d0s=[]
#    omegas12=[]
#    omegas21=[]
#    chis=[]
    natoms1=tnb1.n_atoms
    natoms2=tnb2.n_atoms
    
    resname1=tnb1.top.residue(0).name
    if resname1=='A':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name N6')
    if resname1=='G':
        l_C8_1=tnb1.top.select('name C8')
        l_N6_1=tnb1.top.select('name O6')
    if resname1=='C':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name N4')
    if resname1=='U':
        l_C8_1=tnb1.top.select('name O2')
        l_N6_1=tnb1.top.select('name O4')
    resname2=tnb2.top.residue(0).name
    if resname2=='A':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name N6')
    if resname2=='G':
        l_C8_2=tnb2.top.select('name C8')
        l_N6_2=tnb2.top.select('name O6')
    if resname2=='C':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name N4')
    if resname2=='U':
        l_C8_2=tnb2.top.select('name O2')
        l_N6_2=tnb2.top.select('name O4')

    a2s=tnb2.atom_slice(l_C8_2).xyz
    b2s=tnb2.atom_slice(l_N6_2).xyz
    a1s=tnb1.atom_slice(l_C8_1).xyz
    b1s=tnb1.atom_slice(l_N6_1).xyz
 
    for iframe in range(len(pos1)):
        #compute d0
        CoM1=np.array([sum(pos1[iframe,:,0]),sum(pos1[iframe,:,1]),sum(pos1[iframe,:,2])])/natoms1
        CoM2=np.array([sum(pos2[iframe,:,0]),sum(pos2[iframe,:,1]),sum(pos2[iframe,:,2])])/natoms2
        d0=CoM1-CoM2
        len_d0=np.sqrt(sum(d0**2))
        #compute omega12
        a1=a1s[iframe,0,:]-CoM1
        b1=b1s[iframe,0,:]-CoM1
        c1=np.cross(a1,b1)
        len_c1=np.sqrt(sum(c1**2))
        omega12=np.arccos(np.dot(c1,d0)/(len_c1*len_d0))*180/np.pi
        #compute omega21
        a2=a2s[iframe,0,:]-CoM2
        b2=b2s[iframe,0,:]-CoM2
        c2=np.cross(a2,b2)
        len_c2=np.sqrt(sum(c2**2))
        omega21=np.arccos(np.dot(c2,d0)/(len_c2*len_d0))*180/np.pi
        #compute chi
        chi=np.arccos(np.dot(c1,c2)/(len_c1*len_c2))*180/np.pi

#        len_d0s.append(len_d0)
#        omegas12.append(omegas12)
#        omegas21.append(omegas21)
#        chis.append(chi)
        varis.append(np.array([len_d0,omega12,omega21,chi]))
#    len_d0s=np.array(len_d0s,dtype=np.float32).reshape(len(traj),1)
#    omegas12=np.array(omegas12,dtype=np.float32).reshape(len(traj),1)
#    omegas21=np.array(omegas21,dtype=np.float32).reshape(len(traj),1)
#    chis=np.array(chis,dtype=np.float32).reshape(len(traj),1)
    varis=np.array(varis,dtype=np.float32).reshape(len(traj),4)
#    return np.array([len_d0s,omegas12,omegas21,chis]).T
    return varis
