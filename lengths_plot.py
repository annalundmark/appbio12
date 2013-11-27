#! /usr/bin/env python

'''Script for creating a plot over sequence lengths in a library. 
'''


import argparse
import sys

import collections

import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt

import numpy as np

def main(infile, plotfile, experiment, show_ave): 
	data = file_reader(infile)
	if not data: 
		print "No mapped paired end data"
	else: 
		x, y = length_counter(data)
		average = create_stats(data)
		create_graph(x,y,plotfile, experiment, average, show_ave)
	

def file_reader(infile):
	''' Reads an infile (SAM file) line by line, splits the lines by tabs and 
	stores the fragment length values in a list. '''
	lengths = []
	discarded = 0
	print "Reading file: ", infile
	for ln in open(infile, 'r'): 
		if not ln.startswith("@"): #skips headings
			ln = ln.split('\t')
			if len(ln[9]) != len(ln[10]): 
				print "Not a valid SAM file: length of sequence differs from length of quality"
			if int(ln[8])>0 and check_flag(ln[1]) == True:  
				lengths.append(int(ln[8]))
	print "Done"
	return lengths

def check_flag(flag):
	if bin(int(flag))[-2] == 0: #is the read mapped in a proper pair? 
		return False
	elif len(bin(int(flag))) > 10 and bin(int(flag))[-11] == 1: #does the read pass quality tests? 
		return False
	elif len(bin(int(flag))) > 9 and bin(int(flag))[-10] == 1: #is the read a PCR duplicate? 
		return False
	else:
		return True


def create_stats(lengths): 
	average = int(round(sum(lengths)/float(len(lengths))))
	min_length = min(lengths)
	max_length = max(lengths)
	lower_perc = np.percentile(lengths, 25)
	upper_perc = np.percentile(lengths, 75)
	print "Average length: ", average
	print "Minimum length: ", min_length
	print "Maximum length: ", max_length
	print "25% over ", upper_perc, " bp"
	print "25% under ", lower_perc, " bp"
	return average

def length_counter(lengths): 
	''' Creates two arrays to use when plotting the sequence lengths. 
	One array contains the sequence lengths and one array the number 
	of times that length is observed. '''
	hits_list=[]
	length_list=[]
	hits=collections.Counter(lengths) #creates a dictionary with the sequence lengths as keys and the number of hits as values
	#print hits
	max_length=max(hits.keys())
	#print max_length
	#Create two lists from the hits dictionary, one with ascending sequence lengths and one with corresponding hits
	for i in range (1, max_length+1): 
		x=i
		hits_list.append(x)
		if i in hits: 
			length=hits[i]
			length_list.append(length)
		else: 
			length_list.append(0)
	#print hits_list
	#print length_list
	return hits_list, length_list


def create_graph(x, y, plotfile, experiment, average, show_ave): 
	''' Plots x against y and saves the figure to a file. 
	Also takes in arguments for experiment name as title, average and 
	whether or not to show average length in graph'''
	plt.plot(x, y)
	plt.xlabel("Fragment length (bp)")
	plt.ylabel("Number of fragments")
	if not experiment == None: 
		plt.title(experiment)
	if show_ave: 
		plt.figtext(0.15, 0.05, "Average fragment length: %s bp" % average)
		fig = plt.gcf()
		fig.subplots_adjust(bottom=0.2)
	plt.savefig(plotfile)

	

if __name__ == "__main__": 
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i', '--input', help='Input file (SAM format, required)', required=True)
	parser.add_argument('-p', '--plotfile', help='Name of file to save the plot to (required)', required=True)
	parser.add_argument('-e', '--exp', help='Name of experiment')
	parser.add_argument('-a', '--average', action='store_true', default=False, help='If included, shows average fragment length below graph')
	
	args = parser.parse_args()
	main(args.input, args.plotfile, args.exp, args.average)