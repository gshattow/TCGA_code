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


plot_dir = 'plots/'
OutputFormat = '.png'
parfile = 'parameters.par'
datafile = ''
outfile = ''
ignore_consequence = ''
impact = ''


class ReadData :
	def read_parameters(self, parfile) :
	
		datafile = params.datafile
		outfile = params.outfile
		ignore_consequence = params.ignore_consequence
		impact = params.impact
			
		print datafile, outfile
		return datafile, outfile, ignore_consequence, impact
		
		
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
					if (extra[extra.index('SYMBOL_SOURCE') + 1] == 'HGNC') :
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
		
		print len(genes_all), 'unique genes'
		print len(patients_all), 'unique patients'

		
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
				if item[consequence_index] != ignore_consequence :
					if set(['SYMBOL', 'IND']) <= set(extra) :
						if (extra[extra.index('SYMBOL_SOURCE') + 1] == 'HGNC') & \
							(extra[extra.index('IMPACT') + 1] == impact) :
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

class Analyze :
	def patient_pairs(self, mutation_array) :
		
		print mutation_array.shape
		n_p = (mutation_array.shape)[0]
		n_g = (mutation_array.shape)[1]
		print n_p, n_g
		
		n_pairs = int(n_p/2)
		print n_pairs, 'unique pairs of patients'
		pairs_overlap = np.zeros(n_pairs)
		
		for pp in range(n_pairs) :
			array1 = mutation_array[2*pp]
			array2 = mutation_array[2*pp + 1]
			pair_array = array1 * array2
			pairs_overlap[pp] = np.sum(pair_array)

		print pairs_overlap		

class WriteData :
	def write_mutation_array(self, mutation_array, outfile) :

		np.savetxt(outfile, mutation_array, fmt = '%i')
		
	def write_patients(self, patients_list) :
		np.savetxt('patients.txt', patients_list, fmt = '%s')

	def write_genes(self, genes_list) :
		np.savetxt('genes.txt', genes_list, fmt = '%s')

############################################################



if __name__ == '__main__':
	rd = ReadData()
	az = Analyze()
	wt = WriteData()

	datafile, outfile, ignore_consequence, impact = rd.read_parameters(parfile)
	
	patients_list, genes_list = rd.parse_vep_file(datafile)
	mutation_array = rd.sort_vep_file(datafile, patients_list, genes_list)

	az.patient_pairs(mutation_array)

	wt.write_mutation_array(mutation_array, outfile)
	wt.write_patients(patients_list)
	wt.write_genes(genes_list)