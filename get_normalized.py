import sys
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
import glob
import collections as coll
from collections import Counter
import numpy as np

def get_normalized_pairs(n):
	'''Return a dictionary with keys corresponding to the pairs of residues found 
	within a radius n, and the values to the number of times found in a set of pdb files.\
	This dictionary sets the knowledge of pair-residues at a given frequency found naturally\
	in nature. It is based in 1.110 sequences with known structure with <40% of homology in\
	order to avoid family redundancy. Not necessary for the package.'''
	p = PDBParser(PERMISSIVE=1)
	pdb = glob.glob('./pdbfiles/*.ent')
	pairs = []
	file_list = []	
	
	###### Parsing through PDB files #######
	for filename in pdb:
		s = p.get_structure('X', filename)
		atom_list = np.array([atom for atom in s.get_atoms() if atom.name == 'CB'])
		
			
		if len(atom_list)>2:
			#creates a list containing all atom pairs within a n radius
			ns = Bio.PDB.NeighborSearch(atom_list)
			neighbors = ns.search_all(n)
			file_list.append(filename)
			sys.stderr.write(filename+' processed.\n') #check-point
		else:
			sys.stderr.write(filename+' could not be processed.\n') #check-point
			pass
		
	pairs = [(x.get_parent().get_resname(),y.get_parent().get_resname()) for x,y in neighbors]
	outfile = open( 'normalized_pairs8.py', 'w' )
	counter = dict(Counter(pairs))
	sys.stderr.write(str(len(file_list))+' files processed.\n')			#check-point
	sys.stderr.write('Dictionary length: '+str(len(counter))+'.\n') #check-point
	outfile.write('\nNormalized_pairs_'+str(n)+'='+str(counter))
	outfile.close()



get_normalized_pairs(25)