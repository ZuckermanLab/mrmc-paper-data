#!/usr/bin/env python

import sys
#to install parmed: python setup.py install --home=~
#To use parmed: echo $PYTHONPATH=~/lib/python:$PYTHONPATH
from parmed.amber import *
#import parmed
#recursive. 
#create a dictionary of edges (i,j) such that each edge gives the set of other edges with which it shares an atom
	

#current_path is a list of the bonds on the current recursion, ringset a set containing the indices of bonds in rings 
def find_cycles1(parm, current_path, ring_set_atoms, ring_set_bonds):
	#print "current path: ", current_path
	#need a private copy
	private_current_path = current_path[:]
	iatom = private_current_path[-1]
	if (iatom in ring_set_atoms):
		return ring_set_atoms,ring_set_bonds
	for otheratom in parm.atoms[iatom].bond_partners:
		jatom = otheratom.idx
		#print "looking at bond ",iatom,jatom,private_current_path
		#bond2 can't be the bond we followed to get to bond1
		if ((len(private_current_path)>=2) and (jatom==private_current_path[-2])):
			continue
		if (private_current_path.count(jatom)!=0):
			#we found a ring!
			idx=private_current_path.index(jatom)
			print "found a ring: ",private_current_path[idx:]
			new_ring_set_atoms = ring_set_atoms
			new_ring_set_bonds = ring_set_bonds
			for i in range(idx,len(private_current_path)-1):
				new_ring_set_atoms.add(private_current_path[i])
				new_ring_set_atoms.add(private_current_path[i+1])
				new_ring_set_bonds.add((private_current_path[i],private_current_path[i+1]))
			new_ring_set_atoms.add(iatom)
			new_ring_set_atoms.add(jatom)
			new_ring_set_bonds.add((iatom,jatom))
			return new_ring_set_atoms,new_ring_set_bonds
		else:
			new_list=private_current_path[:]
			new_list.append(jatom)
			new_ring_set_atoms, new_ring_set_bonds = find_cycles1(parm,new_list,ring_set_atoms,ring_set_bonds)
			ring_set_atoms = new_ring_set_atoms
			ring_set_bonds = new_ring_set_bonds
	return ring_set_atoms,ring_set_bonds


def find_cycles(parm,natom):
	ring_set_atoms = set()
	ring_set_bonds = set()
	for iatom in range(0,natom):
		if (iatom not in ring_set_atoms):
			new_ring_set_atoms, new_ring_set_bonds = find_cycles1(parm,[iatom],ring_set_atoms,ring_set_bonds)
			ring_set_atoms = new_ring_set_atoms     
                        ring_set_bonds = new_ring_set_bonds
	return ring_set_atoms,ring_set_bonds


parm = LoadParm(sys.argv[1])
resname = sys.argv[2] #residue name for definitions file
typestart = int(sys.argv[3])
defsfname = sys.argv[4] 
#dihcutoff = float(sys.argv[5])
natom = parm.ptr('natom')
ntypes = parm.ptr('ntypes')
#need to get all the bonds, angles, and dihedrals
nbond = parm.ptr('nbonh') + parm.ptr('mbona')


ring_set_atoms, ring_set_bonds = find_cycles(parm,natom)
print ring_set_atoms
print ring_set_bonds

#to do the definitions file
defsfile = open(defsfname,'w')
defsfile.write('RESI {0}\n'.format(resname))
for iatom in range(0,natom):
	type = iatom + typestart
	defsfile.write('\tATOM {0} {1:d}\n'.format(parm.atoms[iatom].name,type))
for ibond in range(0,nbond):
        iatom = parm.bonds[ibond].atom1         
        jatom = parm.bonds[ibond].atom2
	#determine if the bond is rotatable -- both atoms must be bonded to at least two others and not be part of a ring
	if ((len(iatom.bond_partners)>1) and (len(jatom.bond_partners)>1) and ((iatom.idx,jatom.idx) not in ring_set_bonds) and ((jatom.idx,iatom.idx) not in ring_set_bonds)):
		rotatable = True 
		#changed my mind: every bond not part of the ring should be rotatable (including double bonds)
		#that way we can deal with twisting motions, letting the force field restrain the motion
		#Examine all the dihedrals whose central bond is (iatom,jatom).  If any have k>cutoff, 
		#for dih in iatom.dihedrals:
		#	if (((dih.atom2 is iatom) and (dih.atom3 is jatom)) or ((dih.atom2 is jatom) and (dih.atom3 is iatom))):
		#		if (dih.type.phi_k>dihcutoff):
		#			rotatable = False
		#			break	
	else:
		rotatable = False
	#todo: logic to determine which bonds are rotatable
	if (rotatable):
		defsfile.write('\tBOND {0} {1} ROTATABLE SIDECHAIN\n'.format(iatom.name,jatom.name))
	else:
		defsfile.write('\tBOND {0} {1}\n'.format(iatom.name,jatom.name))
defsfile.write('END\n')
defsfile.close()
