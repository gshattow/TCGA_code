import matplotlib
import scipy
import numpy as np
import pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
import itertools
import time


plot_dir = 'plots/'
OutputFormat = '.png'
datafile = 'Mixed_DNASeq.maf'
# datafile = 'bcm_full.maf'


class ReadData :
	def read_genes(self, datafile) :
		genes_all = []
		patients_all = []


		ref_index = 0
		can1_index = 0
		can2_index = 0
		adj1_index = 0
		adj2_ndex = 0

		for item in file(datafile) :
			item = item.split('\t')
			if len(item) > 1 :
				if item[0] == 'Hugo_Symbol' :
					item = np.array(item)
					chr_index = np.where(item == 'Chromosome')[0][0]
					ref_index = np.where(item == 'Reference_Allele')[0]
					can1_index = np.where(item == 'Tumor_Seq_Allele1')[0]
					can2_index = np.where(item == 'Tumor_Seq_Allele2')[0]
					adj1_index = np.where(item == 'Match_Norm_Seq_Allele1')[0]
					adj2_index = np.where(item == 'Match_Norm_Seq_Allele2')[0]
				if ((item[0] != 'Hugo_Symbol') & (item[9] == 'SNP')):
					genes_all.append(item[0])
					barcode = item[15]
					ids = barcode.split('-')
					patients_all.append(int(ids[2]))
	
		print len(genes_all)
		print len(patients_all)
		print chr_index

		genes_all = list(set(genes_all))
		patients_all = list(set(patients_all))
		
		print len(genes_all)
		print len(patients_all)

		n_g = len(genes_all)
		n_p = len(patients_all)

		cancer = np.zeros((n_p, n_g))
		adjacent = np.zeros((n_p, n_g))
		chromosomes = np.zeros(n_g)

		print cancer.shape

		for item in file(datafile) :
			item = item.split('\t')
			item = np.array(item)
			if len(item) > 1 :
				if ((item[0] != 'Hugo_Symbol') & (item[9] == 'SNP')):
					gene = item[0]
					barcode = item[15]
					ids = barcode.split('-')
					id = int(ids[2])
					cancer_score = 0
					adjacent_score = 0
					if item[ref_index] == item[can1_index] : cancer_score = cancer_score + 1
					if item[ref_index] == item[can2_index] : cancer_score = cancer_score + 1
					if item[ref_index] == item[adj1_index] : adjecent_score = adjacent_score + 1
					if item[ref_index] == item[adj2_index] : adjecent_score = adjacent_score + 1
					p_index = patients_all.index(id)
					g_index = genes_all.index(gene)
					cancer[p_index, g_index] = cancer_score
					adjacent[p_index, g_index] = adjacent_score
					if ((item[chr_index] == 'X') | (item[chr_index] == 'Y')) :
						chromosomes[g_index] = 23
					else : chromosomes[g_index] = item[chr_index]
	
		ii = np.argsort(chromosomes)
		
		chromosomes = chromosomes[ii]
		genes_all = np.array(genes_all)
		genes_all = list(genes_all[ii])
		cancer = cancer[:, ii]
		adjacent = adjacent[:, ii]

	
		outfile = 'cancer.txt'
		np.savetxt(outfile, cancer, delimiter='\t', fmt = '%d')

		genefile = 'genes.txt'
		np.savetxt(genefile, genes_all, delimiter='\n', fmt = '%s')

		patientfile = 'patients.txt'
		np.savetxt(patientfile, patients_all, delimiter='\n', fmt = '%s')


		G1INDEX = genes_all.index('PBRM1')
		print 'PBRM1 gene index:', G1INDEX
		G2INDEX = genes_all.index('ARID1A')
		G3INDEX = genes_all.index('SMARCA4')

		print np.sum(cancer[:, G1INDEX]), 'mutations in the PBRM1 gene'
		print np.sum(cancer[:, G2INDEX]), 'mutations in the ARID1A gene'
		print np.sum(cancer[:, G3INDEX]), 'mutations in the SMARCA4 gene'

		return  patients_all, genes_all, cancer, adjacent, chromosomes

class Results :
	def Cross_Correlation(self, N, cancer) :


		# populate an array with gene pairs for each patient
		# i.e. if a patient has mutations in genes 1 and 3
		# we add 1 to network[1,3]

		tic = time.clock()
		network = np.zeros((N,N))
		print len(cancer)
		for pp in range(len(cancer)) :	
			ww = np.where(cancer[pp] > 0)[0]	
			for pair in itertools.combinations(ww, 2):
				network[pair[0], pair[1]] = network[pair[0], pair[1]] + 1

		toc = time.clock()
		print toc - tic
		
		print network.shape
		return network

class PlotData :
	def N_mutations_per_patient(self, n_p, n_g, cancer_mutation_count, adjacent_mutation_count) :
		fig = plt.figure(figsize=(10, 3))
		ax = fig.add_subplot(111)

		index = np.arange(n_p)
		bar_width = 0.25

		opacity = 0.4

		cancer_bars = plt.bar(index, cancer_mutation_count, bar_width,
						 alpha=opacity,
						 color='b',
						 label='Cancer')

		adjacent_bars = plt.bar(index, adjacent_mutation_count, bar_width,
						 alpha=opacity,
						 color='r',
						 label='Adjacent')


		plt.xlabel('Patients')
		plt.ylabel('Mutation')
		plt.xticks(index + bar_width)
		ax.set_xticklabels(patients_all, rotation=45 ) ;

		plt.legend()

		outputFile = plot_dir + 'N_mutations_per_patient' + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()

	def N_patients_with_N_mutations(self, cancer_mutation_count) :
		fig = plt.figure(figsize=(7, 6))
		ax = fig.add_subplot(111)

		bins = np.max(cancer_mutation_count) - np.min(cancer_mutation_count)
		hist, bin_edges = np.histogram(cancer_mutation_count, bins)
		bin_centres = bin_edges + 0.5
		bin_centres = bin_centres[0:bins]

		ax.plot(bin_centres, hist)

		plt.xlabel('Number of Mutations')
		plt.ylabel('Number of Patients')

		outputFile = plot_dir + 'N_mutations_histogram' + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()

	def Network(self, min_pairs, network, mutations_per_gene, chromosomes) :
	
		tic1 = time.clock()
		fig = plt.figure(figsize=(10, 9))
		ax = fig.add_subplot(111)
		plt.axis([-1.1,1.1,-1.1,1.1])
		ax.set_xticklabels([])
		ax.set_yticklabels([])

		ww = np.where(network > min_pairs-1)
#		print ww
		N_pairs = len(ww[0])
		print N_pairs, 'pairs occur at least', min_pairs, 'times'
#		print network[np.where(network > min_pairs-1)]
		
		uniq_genes = list(set(np.ravel(ww)))
		print len(uniq_genes), 'genes are involved in the pairs'
		N_genes = len(uniq_genes)
		
		theta = 2*np.pi*np.linspace(0, 1, N_genes)
		r = 1. 
		x = r*np.sin(theta)
		y = r*np.cos(theta)
		colors = chromosomes[uniq_genes]/23.
		area = mutations_per_gene[uniq_genes]
		print 'max area:', np.max(area)


		print area
	
		tic = time.clock()
		plt.scatter(x, y, s=area, c=colors, alpha=0.5)
		for gg in range(N_pairs):
			g1 = uniq_genes.index(ww[0][gg])
			g2 = uniq_genes.index(ww[1][gg])
#			print gg, g1, g2, ww[0][gg], ww[1][gg], network[ww[0][gg], ww[1][gg]]
			ax.plot((x[g1], x[g2]), (y[g1], y[g2]), c = 'k', \
				alpha = network[ww[0][gg], ww[1][gg]]/np.max(network))

		outputFile = plot_dir + 'gene_network_' + str(N_pairs) + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()
		toc = time.clock()
		print toc - tic
		print 'total time to make figure:', toc - tic1

	def Network_by_chromosome(self, cancer, chromosomes) :

		fig = plt.figure(figsize=(10, 9))
		ax = fig.add_subplot(111)
		plt.axis([-1.1,1.1,-1.1,1.1])
		ax.set_xticklabels([])
		ax.set_yticklabels([])

		N = 100
		
		for ii in range(23) : print len(np.where(chromosomes == ii + 1)[0])
		test = cancer[:, chromosomes == 1]
		print test.shape
		test[test > 1] = 1
		
		
		
if __name__ == '__main__':
	rd = ReadData()
	res = Results()
	pd = PlotData()

	patients_all, genes_all, cancer, adjacent, chromosomes = rd.read_genes(datafile)
	
	cancer_mutation_count = np.sum(cancer, axis = 1)
	mutations_per_gene = np.sum(cancer, axis = 0)

	adjacent_mutation_count = np.sum(adjacent, axis = 1)
	# print adjacent_mutation_count
	print sum(adjacent_mutation_count)

	n_p = len(patients_all)
	n_g = len(genes_all)
	
	N = n_g
	network = res.Cross_Correlation(N, cancer)
	
#	pd.N_mutations_per_patient(n_p, n_g, cancer_mutation_count, adjacent_mutation_count)
#	pd.N_patients_with_N_mutations(cancer_mutation_count)

	min_pairs = 4
	pd.Network(min_pairs, network, mutations_per_gene, chromosomes)
#	pd.Network_by_chromosome(cancer, chromosomes)