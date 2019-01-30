#!/usr/bin/env python

import sys
#to install parmed: python setup.py install --home=~
#To use parmed: echo $PYTHONPATH=~/lib/python:$PYTHONPATH
from parmed.amber import *
import networkx as nx

#import parmed
#recursive. 
#create a dictionary of edges (i,j) such that each edge gives the set of other edges with which it shares an atom
	

parm = LoadParm(sys.argv[1])
resname = sys.argv[2] #residue name for definitions file
typestart = int(sys.argv[3])
defsfname = sys.argv[4] 
#dihcutoff = float(sys.argv[5])
natom = parm.ptr('natom')
ntypes = parm.ptr('ntypes')
#need to get all the bonds, angles, and dihedrals
nbond = parm.ptr('nbonh') + parm.ptr('mbona')

molgraph=nx.Graph()
molgraph.add_nodes_from([0,natom])
for ibond in range(0,nbond):
        iatom = parm.bonds[ibond].atom1.idx
        jatom = parm.bonds[ibond].atom2.idx
	molgraph.add_edge(iatom,jatom)

cycles= nx.cycle_basis(molgraph)
ring_set_atoms=set()
ring_set_bonds=set()
print cycles
for cycle in cycles:
	n=len(cycle)
	for i in range(0,n):
		iatom=cycle[i]
		if (i==n-1):
			jatom=cycle[0]
		else:
			jatom=cycle[i+1]
		ring_set_atoms.add(iatom)
		ring_set_bonds.add((iatom,jatom))


#ring_set_atoms, ring_set_bonds = find_cycles(parm,natom)
print ring_set_atoms
print ring_set_bonds

#sys.exit(0)
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
