###############################################################################
### functions to identify the indexes of
### the atoms to define the backbone dihedrals
### in an RNA molecule
###############################################################################
### all the definitions are taken from 
### from http://x3dna.org/highlights/torsion-angles-of-nucleic-acid-structures
###############################################################################
### written by Giovanni Pinamonti, PhD student @ SISSA, Trieste, Italy
### date: November 2015
###############################################################################

import mdtraj

def reference_dih():
    bb_list=["O3'","P","O5'","C5'","C4'","C3'"]
    bb_circ_list=bb_list+bb_list[:3]
    res_ndx_list=[-1,0,0,0,0,0,0,+1,+1]
    ### print the definitions of dihedrals
    name_list=['alpha','beta','gamma','delta','espilon','zeta']
    print 'Backbone dihedrals definitions:'
    for i in range(len(bb_circ_list)-3):
        print name_list[i],zip(bb_circ_list,res_ndx_list)[i:i+4]
    print ''
    sr_list=["C4'","O4'","C1'","C2'","C3'"]
    sr_circ_list=sr_list+sr_list[:3]
    name_list=['v0','v1','v2','v3','v4']
    print 'Sugar ring torsional angle'
    for i in range(len(sr_circ_list)-3):
        print name_list[i],sr_circ_list[i:i+4]
    print ''

    R_list=['A','G']
    Y_list=['U','C']
    Y_atoms=["O4'","C1'","N1","C2"]
    R_atoms=["O4'","C1'","N9","C4"]
    print 'Chi torsional angles:'
    print "for pyrimidines(Y): O4'-C1'-N1-C2"
    print "for purines(R)    : O4'-C1'-N9-C4"


    print '(from http://x3dna.org/highlights/torsion-angles-of-nucleic-acid-structures )'



# Find index of backbone atoms 
#      alpha:   O3'(i-1)-P-O5'-C5'
#      beta:    P-O5'-C5'-C4'
#      gamma:   O5'-C5'-C4'-C3'
#      delta:   C5'-C4'-C3'-O3'
#      epsilon: C4'-C3'-O3'-P(i+1)
#      zeta:    C3'-O3'-P(i+1)-O5'(i+1)
# (from http://x3dna.org/highlights/torsion-angles-of-nucleic-acid-structures)
def get_dihedrals_ndx(top):
    bb_list=["O3'","P","O5'","C5'","C4'","C3'"]
    bb_circ_list=bb_list+bb_list[:3]
    res_ndx_list=[-1,0,0,0,0,0,0,+1,+1]
    name_list=['alpha','beta','gamma','delta','espilon','zeta']
    # index search and store
    dihedrals_list=[]
    for ires in range(top.n_residues):
        tmp_dih_list=[]
        for i in range(len(bb_circ_list)-3):
            tmp_list=[]
            for (name,i_ndx) in zip(bb_circ_list,res_ndx_list)[i:i+4]:
                if ires+i_ndx>top.n_residues-1:
                    continue
                for atom in top.residue(ires+i_ndx).atoms_by_name(name):
                    tmp_list.append(atom.index)
            tmp_dih_list.append(tmp_list)
        dihedrals_list.append(tmp_dih_list)
    
    list_of_dihedrals=[[],[],[],[],[],[]]
    for i in range(top.n_residues):
        iname=0
        for (name,ndx) in zip(name_list,dihedrals_list[i]):
            if len(ndx)==4:
                list_of_dihedrals[iname].append(ndx)
            iname+=1

    return list_of_dihedrals


# Sugar conformational parameters: 
#      v0: C4'-O4'-C1'-C2'
#      v1: O4'-C1'-C2'-C3'
#      v2: C1'-C2'-C3'-C4'
#      v3: C2'-C3'-C4'-O4'
#      v4: C3'-C4'-O4'-C1'
# (from http://x3dna.org/highlights/torsion-angles-of-nucleic-acid-structures)
def get_pucker_ndx(top):
    sr_list=["C4'","O4'","C1'","C2'","C3'"]
    sr_circ_list=sr_list+sr_list[:3]
    name_list=['v0','v1','v2','v3','v4']
    pucker_list=[]
    for ires in range(top.n_residues):
        tmp_pck_list=[]
        for i in range(len(sr_list)):
            tmp_list=[]
            for name in sr_circ_list[i:i+4]:
                for atom in top.residue(ires).atoms_by_name(name):
                    tmp_list.append(atom.index)
            tmp_pck_list.append(tmp_list)
        pucker_list.append(tmp_pck_list)
    
    list_of_pucker=[[],[],[],[],[]]
    for i in range(top.n_residues):
        iname=0
        for (name,ndx) in zip(name_list,pucker_list[i]):
            if len(ndx)==4:
                list_of_pucker[iname].append(ndx)
            iname+=1
 
    return list_of_pucker

# Chi torsional angles:
# for pyrimidines(Y): O4'-C1'-N1-C2
# for purines(R)    : O4'-C1'-N9-C4
# (from http://x3dna.org/highlights/torsion-angles-of-nucleic-acid-structures)
def get_chi_ndx(top):
    R_list=['A','G']
    Y_list=['U','C']
    Y_atoms=["O4'","C1'","N1","C2"]
    R_atoms=["O4'","C1'","N9","C4"]
    chi_list=[]
    for ires in range(top.n_residues):
        res=top.residue(ires)
        tmp_chi_list=[]
        lista=[]
        if res.name in R_list:
            lista=R_atoms
        if res.name in Y_list:
            lista=Y_atoms
        for name in lista:
            for atom in res.atoms_by_name(name):
                tmp_chi_list.append(atom.index)
        if len(tmp_chi_list)==4:
            chi_list.append(tmp_chi_list)
    return chi_list
