#!/usr/bin/env python3

import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
import glob
import collections as coll
from collections import Counter
import numpy as np
from .normalized_pairs import *

def get_info(filename):
	'''
	Return header. Function adapted from Biopython Package.\n
	get_info(filename)\n
	Filename needs to be a PDB file format (*.ent or *.pdb)
	'''
	p = PDBParser(QUIET=True)
	s = p.get_header()

def parse_atoms_infile(filename):
	'''
	Parse a PDB file and return atom list.\n
	parse_atoms_infile(filename):\n
	File needs to be a PDB file format (*.ent or *.pdb)
	'''
	p = PDBParser(QUIET=True)
	s = p.get_structure("X", filename)
	atom_list = [atom for atom in s.get_atoms() if atom.name == 'CB']
	return atom_list

def get_neighbors(filename,r):
	'''
	Return atom pairs with a neighbor-joining algorithm at a given radius 'r'.\n
	get_neighbors(filename,r)
	'''
	atom_list = parse_atoms_infile(filename)

	if len(atom_list)>2:
		dictio = {}
		ns = Bio.PDB.NeighborSearch(atom_list)
		ns_2= ns.search_all(r)
		return ns_2
	else:
		ns_2 = []
		return ns_2


def get_position(filename,r):
	'''
	Return position list of the residue-pairs at a given radius 'r'.\n
	get_position(filename,r)
	'''
	ns_2 = get_neighbors(filename,r)
	pos1 = []
	pos = []
	residues = [list(x.get_full_id()) for x,y in ns_2]
	for element in list(residues):
		pos1.append(list(element[3]))
		
	for element in pos1:
		pos.append(element[1])
	return pos

	
def get_infile_pairs(filename,r):
	'''
	Return dictionary key{residue pair}:value{number of occurrences in inFile}
	at a given radius 'r'.\n
	get_infile_pairs(filename,r)
	'''
	ns_2 = get_neighbors(filename, r)
	pos = get_position(filename, r)
	reduced = []
	counter = {}
	
	if len(ns_2)>2:
		pairs = [(x.get_parent().get_resname(),y.get_parent().get_resname()) for x,y in ns_2]
		pairs_pos = list(zip(pairs,pos))
		for element in list(pairs_pos):
			reduced.append(element[1])
		
		number = dict(Counter(reduced))
		
		for element in pairs_pos:
			for key in number:
				if element[1] == key:
					counter[element] = number.get(key)
		return counter

	else:
		pass
						

def compare_normalized_pairs(filename, r):
	'''
	Return dictionary key{residue pair}:value{number of occurrences 
	in data normalization at a given radius 'r'}. 'r' should be 
	chosen between 5, 10, 15, 20 and 25Å in order to be comparable 
	with the normalized dictionaries.\n
	compare_normalized_pairs(filename, r)
	'''
	dictio_normalized = {}
	inp = get_infile_pairs(filename, r)
	for k,v in inp:
		if r == 5:
			radius = Normalized_pairs_5
			if k in radius.keys():
				dictio_normalized[k] = radius.setdefault(k)			
		elif r == 15:
			radius = Normalized_pairs_15
			if k in radius.keys():
				dictio_normalized[k] = radius.setdefault(k)
		elif r == 20: 
			radius = Normalized_pairs_20
			if k in radius.keys():
				dictio_normalized[k] = radius.setdefault(k)
		elif r == 25:
			radius = Normalized_pairs_25
			if k in radius.keys():
				dictio_normalized[k] = radius.setdefault(k)
		else:
			radius = Normalized_pairs_10
			if k in radius.keys():
				dictio_normalized[k] = radius.setdefault(k)
					
	return dictio_normalized

def dictio_pairs(filename,r):
	'''Return dictionary 
	key{residue pair tuple, and position}:value{statistical potential energy}
	at a given radius 'r'. 'r' should be chosen between 5, 10, 15, 20 and 25Å in order
	to be comparable with the normalized dictionaries.\n
	dictio_pairs(filename, r)
	'''
	kb = 0.0019872041  #Boltzmann constant
	T = 297.15  #Temperature
	ct = np.multiply(kb,T)
	dictio_infile = get_infile_pairs(filename, r)
	dictio_normalized = compare_normalized_pairs(filename, r)
	dictio_pairs = {}
	list_tupla = []
	norm = dictio_normalized.keys()
	
	for key1,key2 in dictio_infile:
		tupla = (key1,key2)
		if key1 in dictio_normalized.keys():
			quocient = np.divide(dictio_infile.get(tupla),dictio_normalized.get(key1))
			ln = np.log(quocient)
			U = float(-(np.multiply(ct,ln)))
			dictio_pairs[tupla] = U
			
	return dictio_pairs


def get_potential(filename,r):
	'''
	Return dictionary key{position-residue}:value{statistical potential per position}
	at a given radius 'r'. 'r' should be chosen between 5, 10, 15, 20 and 25Å in order
	to be comparable with the normalized dictionaries.\n
	get_potential(filename, r)
	'''
	dictio = dictio_pairs(filename, r)
	positions_set = set()
	lista_todo = []
	dictio_pos_res = {}
	dictio_freqs = {}
	residue = []
	final = []
	list_final = []
	sum = []
		
	for k,v in dictio:
		tupla = (k,v)
		positions_set.add(v)				# set with the different positions
		element = [k,v,dictio.get(tupla)]
		lista_todo.append(element)		#list of lists (pair, position, frequency)

	positions_list = list(positions_set)
	lista_todo.sort(key=lambda x: x[1])	
	
	for element in lista_todo:
		dictio_pos_res[element[1]] = str(element[1])+'-'+str(list(element[0])[0])
		# dictionary that contains key = position, and value = 1st residue + position
	
	for element in lista_todo:
		tupletix = (element[1],element[2])
		residue.append(tupletix)

	for position in positions_list:
		for element in residue:
			if position == element[0]:
				final.append(element[1])
		list_final.append(final)
		final = []
	for element in list_final:
		final = float(np.sum(element))
		sum.append(final)
	
	dictio_freqs = dict(zip(positions_list, sum))
	
	return dictio_freqs
