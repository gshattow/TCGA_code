import matplotlib
import scipy
import numpy as np
import pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

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

		genes_all = list(set(genes_all))
		patients_all = list(set(patients_all))

		print len(genes_all)
		print len(patients_all)

		n_g = len(genes_all)
		n_p = len(patients_all)

		cancer = np.zeros((n_p, n_g))
		adjacent = np.zeros((n_p, n_g))

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
	
		outfile = 'cancer.txt'
		np.savetxt(outfile, cancer, delimiter='\t', fmt = '%d')

		genefile = 'genes.txt'
		np.savetxt(genefile, genes_all, delimiter='\n', fmt = '%s')

		patientfile = 'patients.txt'
		np.savetxt(patientfile, patients_all, delimiter='\n', fmt = '%s')



		return  patients_all, genes_all, cancer, adjacent




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
		print bin_edges
		bin_centres = bin_edges + 0.5
		bin_centres = bin_centres[0:bins]
		print bin_centres

		ax.plot(bin_centres, hist)

		plt.xlabel('Number of Mutations')
		plt.ylabel('Number of Patients')

		outputFile = plot_dir + 'N_mutations_histogram' + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()


if __name__ == '__main__':
	rd = ReadData()
	pd = PlotData()

	patients_all, genes_all, cancer, adjacent = rd.read_genes(datafile)
	
	cancer_mutation_count = np.sum(cancer, axis = 1)
	print cancer_mutation_count

	adjacent_mutation_count = np.sum(adjacent, axis = 1)
	# print adjacent_mutation_count
	print sum(adjacent_mutation_count)

	n_p = len(patients_all)
	n_g = len(patients_all)
	
	pd.N_mutations_per_patient(n_p, n_g, cancer_mutation_count, adjacent_mutation_count)
	pd.N_patients_with_N_mutations(cancer_mutation_count)
