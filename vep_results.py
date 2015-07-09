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

SIFT = "%03d" % int(params.sift*100)


class ReadData :		
	def read_mutation_array(self) :
		file = 'mutation_array_c' + params.chrom + '_s' + SIFT + '.dat'
		mutation_array = np.loadtxt(file, unpack=True)
		mutation_array = np.transpose(mutation_array)
		
		return mutation_array
		

	def read_patients(self) :
		file = 'patients_c' + params.chrom + '.txt'
		patients_list = np.loadtxt(file, dtype=np.str)
		
		return patients_list

	def read_genes(self) :
		file = 'genes_c' + params.chrom + '.txt'
		genes_list = np.loadtxt(file, dtype=np.str)
		
		return genes_list
			

class Results :
	def patient_pairs(self, mutation_array) :
		
		tic = time.clock()
		print mutation_array.shape
		n_p = (mutation_array.shape)[0]
		n_g = (mutation_array.shape)[1]
		print n_p, 'people', n_g, 'genes'
		
		b_mutation_array = np.where(mutation_array > 0, 1., 0,)
		print b_mutation_array
		n_pairs = int(n_p/2)
		print n_pairs, 'unique pairs of patients'
		n_runs = params.n_runs
		pairs_overlap = np.zeros((n_runs, n_pairs))

		
		list_p = np.linspace(0, n_p - 1, n_p)
		for run in range(n_runs) :
			randomized_list = sample(list_p, n_p)
			for pq in range(n_pairs) :
				array1 = b_mutation_array[randomized_list[2*pq]]
				array2 = b_mutation_array[randomized_list[2*pq + 1]]
				pair_array = array1 * array2
				pairs_overlap[run][pq] = np.sum(pair_array)

#		print 'overlapping pairs', pairs_overlap		
		toc = time.clock()
		print toc - tic, 'seconds to pair patients'
		
		return pairs_overlap

	def gene_pairs(self, mutation_array) :

		n_p = (mutation_array.shape)[0]
		n_g = (mutation_array.shape)[1]

		# populate an array with gene pairs for each patient
		# i.e. if a patient has mutations in genes 1 and 3
		# we add 1 to network[1,3]

		tic = time.clock()
		gene_pair_array = np.zeros((n_g,n_g))
		print len(mutation_array)
		for pp in range(n_p) :	
			ww = np.where(mutation_array[pp] > 0)[0]	
			for pair in itertools.combinations(ww, 2):
				gene_pair_array[pair[0], pair[1]] += 1

		toc = time.clock()
		print toc - tic, 'seconds to cross correlate pairs'
		
		print gene_pair_array.shape
		
		print gene_pair_array
		return gene_pair_array


class WriteData :
		
	def write_pairs_overlap(self, pairs_overlap) :
		file = 'patient_pairs_overlap_c' + params.chrom + '_s' + SIFT + '.txt'
		np.savetxt(file, pairs_overlap, fmt = '%i')

	def write_gene_pairs(self, gene_pair_array) :
		file = 'gene_pairs_c' + params.chrom + '_s' + SIFT + '.txt'
		np.savetxt(file, gene_pair_array, fmt = '%i')

############################################################



if __name__ == '__main__':
	rd = ReadData()
	res = Results()
	wt = WriteData()

	
	mutation_array = rd.read_mutation_array()
	patients_list = rd.read_patients()
	genes_list = rd.read_genes()
	
	pairs_overlap = res.patient_pairs(mutation_array)
	wt.write_pairs_overlap(pairs_overlap)

	gene_pair_array = res.gene_pairs(mutation_array)
	wt.write_gene_pairs(gene_pair_array)