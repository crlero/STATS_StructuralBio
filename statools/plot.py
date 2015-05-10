#!/usr/bin/env python3

import sys
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from pylab import *
from . import functions

def make_plot(filename,r=10):
	'''
	Show, and save plot of statistical potentials
	make_plot(filename, r) -> where filename needs to be a PDB file format
	(*.pdb or the old version *.ent), and 'r' a given radius (5, 10, 15, 20, or
	25Å, for the neighbor-joining algorithm)
	'''
	
	# data sets
	dictio = functions.get_potential(filename, r)
	y = []
	residues = sorted(dictio)	# residues
	for element in residues:
		y.append(dictio.get(element))	#potential energy
	x = list(range(len(residues)))
	length = len(residues)
	
	# plotting
	fig = plt.figure(length,edgecolor='k',figsize=(4*3.13,2*3.13))
	gs = GridSpec(1,1,bottom=0.20,left=0.07,right=0.93,top=0.90)
	ax = fig.add_subplot(111)
	fig.subplots_adjust(bottom=0.2)
	fig.patch.set_visible(False)
	ax.plot(x,y, 'bo', x, y, 'k', markerfacecolor='m', color='gray')
	fig.patch.set_visible(False)
	
	# draw a default hline at y=1 that spans the xrange
	l = ax.axhline(y=0,color='r')
	# draw a default vline at x=1 that spans the yrange
	l = ax.axvline(x=0, color='r')

	# definitions for the axes
	plt.xticks(x,residues)
	xlabels = residues
	# Changing the label rotation from the x axis
	ax.set_xticklabels(xlabels,size=6,rotation=45) 
	plt.xlim(0,length)
	
	# axis labels customization
	ttext = plt.title("Statistical Potential for '"+filename+"'")
	ytext = plt.ylabel('Potential Energy (kcal/mol)')
	xtext = plt.xlabel('Position-Residues')
	plt.setp(ttext, size='large', color='r', weight='bold', style='normal')
	plt.setp(xtext, size='medium', name='Bitstream Vera Sans', weight='light')
	plt.setp(ytext, size='medium', name='Bitstream Vera Sans', weight='light')
	
	# save image as png
	plt.savefig(filename+'.png', dpi=200)
	
	#show stderr
	plt.show()

