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


plot_dir = 'plots/'
OutputFormat = '.png'
datafile = ''
outfile = ''
ignore_consequence = ''
impact = ''


class ReadData :		
		
	def parse_vep_file(self, datafile) :

		genes_all = []
		patients_all = []

		location_index = 0
		extra_index = 0
		gene_index = 0
		patient_index = 0
		consequence_index = 0
		genes_all = []
		patients_all = []
		chromosomes = []
		
		ii = 0
		tic = time.clock()
		for item in file(datafile) :
			ii = ii + 1
			if item.startswith('##'):
				continue
			elif item.startswith('#') :
				print 'line number', ii, 'has the index names'
				item = item.split()
				consequence_index = item.index('Consequence')
				location_index = item.index('Location')
				extra_index = item.index('Extra')
			else :
				item = item.split()
				loc = item[location_index].split(':')
				extra = item[extra_index].replace(';',' ').replace('=', ' ').split()
				if set(['SYMBOL', 'IND']) <= set(extra) :
					if (extra[extra.index('SYMBOL_SOURCE') + 1] == params.catalogue) :
						gene_index = extra.index('SYMBOL') + 1
						patient_index = extra.index('IND')  + 1
						genes_all.append(extra[gene_index])
						patients_all.append(extra[patient_index])

		toc = time.clock()	
		print toc - tic, 'seconds to read in file, ', \
			float(len(genes_all))/(toc-tic), 'per second'		
					
		print len(genes_all), 'genes read'
		print len(patients_all), 'patients read'

		genes_all = list(set(genes_all))
		patients_all = list(set(patients_all))
		
		print len(genes_all), 'intragenic mutations read'

		
#		print genes_all
#		print patients_all
		
		return patients_all, genes_all
		
	def sort_vep_file(self, datafile, patients_list, genes_list) :
		
		n_g = len(genes_list)
		n_p = len(patients_list)
		mutation_array = np.zeros((n_p, n_g))
		
		tic = time.clock()
		for item in file(datafile) :
			if item.startswith('##'):
				continue
			elif item.startswith('#') :
				item = item.split()
				consequence_index = item.index('Consequence')
				location_index = item.index('Location')
				extra_index = item.index('Extra')
			else :
				item = item.split()
				loc = item[location_index].split(':')
				extra = item[extra_index].replace(';',' ').replace('=', ' ').split()
				if item[consequence_index] != params.ignore_consequence :
					if set(['SYMBOL', 'IND', 'SIFT']) <= set(extra) :
						if (extra[extra.index('SYMBOL_SOURCE') + 1] == params.catalogue) :
#							(extra[extra.index('IMPACT') + 1] == params.impact) :
							sift = extra[extra.index('SIFT') + 1].replace('(', ' ').replace(')', ' ').split()
							if float(sift[1]) < 0.05 :
								gene_index = extra.index('SYMBOL') + 1
								patient_index = extra.index('IND')  + 1
								ii_p = patients_list.index(extra[patient_index])
								ii_g = genes_list.index(extra[gene_index])
								mutation_array[ii_p][ii_g] += 1

		toc = time.clock()	
		print toc - tic, 'seconds to read in file, ', \
			float(len(genes_list))/(toc-tic), 'per second'		
			
#		print mutation_array
		
		return mutation_array


class WriteData :
	def write_mutation_array(self, mutation_array) :
		file = 'mutation_array_' + params.chrom + '_' + params.impact + '.dat'
		np.savetxt(file, mutation_array, fmt = '%i')
		
	def write_patients(self, patients_list) :
		file = 'patients_' + params.chrom + '_' + params.impact + '.txt'
		np.savetxt(file, patients_list, fmt = '%s')

	def write_genes(self, genes_list) :
		file = 'genes_' + params.chrom + '_' + params.impact + '.txt'
		np.savetxt(file, genes_list, fmt = '%s')
		

############################################################



if __name__ == '__main__':
	rd = ReadData()
	res = Results()
	wt = WriteData()

	
	datafile = params.datafile + params.chrom
	patients_list, genes_list = rd.parse_vep_file(datafile)
	mutation_array = rd.sort_vep_file(datafile, patients_list, genes_list)


	wt.write_mutation_array(mutation_array)
	wt.write_patients(patients_list)
	wt.write_genes(genes_list)
