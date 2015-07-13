import matplotlib
import scipy
import numpy as np
import pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
import itertools
import time
import params
from random import sample
import os
import subprocess


plot_dir = 'plots/'
OutputFormat = '.png'

SIFT = "%03d" % int(params.sift*100)


class ReadData :			
	
	def parse_vep_file(self, datafile) :

		genes_all = []
		individuals_all = []

		location_index = 0
		extra_index = 0
		gene_index = 0
		individual_index = 0
		consequence_index = 0
		genes_all = []
		individuals_all = []
		chromosomes = []
		
		
		ii = 0
		tic = time.clock()
		for item in file(datafile) :
			ii = ii + 1
			if item.startswith('#') :
				print 'line number', ii, 'has the index names'
				item = item.split()
				location_index = item.index('Location')
				extra_index = item.index('Extra')
			else :
				item = item.split()
				loc = item[location_index].split(':')
				extra = item[extra_index].replace(';',' ').replace('=', ' ').split()
				if (extra[extra.index('SYMBOL_SOURCE') + 1] == params.catalogue) :
					gene_index = extra.index('SYMBOL') + 1
					individual_index = extra.index('IND')  + 1
					genes_all.append(extra[gene_index])
					individuals_all.append(extra[individual_index])

		toc = time.clock()	
		print toc - tic, 'seconds to read in file, ', \
			float(len(genes_all))/(toc-tic), 'per second'		
					
		print len(genes_all), 'intragenic mutations read'

		genes_all = list(set(genes_all))
		individuals_all = list(set(individuals_all))
		
		print len(genes_all), 'unique genes'
		print len(individuals_all), 'unique individuals'

				
		return individuals_all, genes_all, location_index, extra_index
		
	def sort_vep_file(self, datafile, individuals_list, genes_list, location_index, extra_index) :
		
		n_g = len(genes_list)
		n_p = len(individuals_list)
		mutation_array = np.zeros((n_p, n_g))
		
		tic = time.clock()
		for item in file(datafile) :
			if item.startswith('r') :
				item = item.split()
				loc = item[location_index].split(':')
				extra = item[extra_index].replace(';',' ').replace('=', ' ').split()
				if (extra[extra.index('SYMBOL_SOURCE') + 1] == params.catalogue) :
					sift = extra[extra.index('SIFT') + 1].replace('(', ' ').replace(')', ' ').split()
					if float(sift[1]) < params.sift :
						gene_index = extra.index('SYMBOL') + 1
						individual_index = extra.index('IND')  + 1
						ii_p = individuals_list.index(extra[individual_index])
						ii_g = genes_list.index(extra[gene_index])
						mutation_array[ii_p][ii_g] += 1

		toc = time.clock()	
		print toc - tic, 'seconds to read in file, ', \
			float(len(genes_list))/(toc-tic), 'per second'		
			
		
		return mutation_array


class WriteData :
	def write_mutation_array(self, mutation_array) :
		file = 'mutation_array_c' + params.chrom + '_s' + SIFT + '.dat'
		np.savetxt(file, mutation_array, fmt = '%i')
		
	def write_individuals(self, individuals_list) :
		file = 'individuals_c' + params.chrom + '.txt'
		np.savetxt(file, individuals_list, fmt = '%s')

	def write_genes(self, genes_list) :
		file = 'genes_c' + params.chrom + '.txt'
		np.savetxt(file, genes_list, fmt = '%s')
		

############################################################



if __name__ == '__main__':
	rd = ReadData()
	wt = WriteData()

	
	datafile = params.datafile
		
	individuals_list, genes_list, location_index, extra_index = \
		rd.parse_vep_file(datafile)
	mutation_array = \
		rd.sort_vep_file(datafile, individuals_list, genes_list, location_index, extra_index)


	wt.write_mutation_array(mutation_array)
	wt.write_individuals(individuals_list)
	wt.write_genes(genes_list)
