import numpy as np
import scipy
import matplotlib
import params
import pylab
import matplotlib.pyplot as plt
from random import sample
import matplotlib.colors as colors
import matplotlib.cm as cm

font = {'family':'serif',
	'size':20	}
plt.rc('font', **font)


OutputFormat = '.png'
chrom = params.chrom
impact = params.impact


class PlotData :		

	def patient_overlap(self) :
		pairs_overlap = []
		patient_pair_file = 'patient_pairs_overlap_' + chrom + '_' + impact + '.txt'
		for item in file(patient_pair_file) :
			item = item.split()
			pairs_overlap.append(map(float, item))
		
		
		pairs_overlap = np.array(pairs_overlap)		
		print pairs_overlap
		print pairs_overlap.shape

		min_p = np.min(pairs_overlap)
		max_p = np.max(pairs_overlap)
		print 'min, max = ', min_p, max_p
		nbins = int(max_p) + 1
		n_runs = len(pairs_overlap)


		print np.min(pairs_overlap), np.max(pairs_overlap)
		nbins = int(np.max(pairs_overlap))
		bin_centres = np.linspace(0, nbins, nbins)
		bin_edges = np.linspace(-0.5, nbins + 0.5, nbins + 1)


		fig = plt.figure(frameon=False, figsize=(10, 9))
		ax = fig.add_subplot(111)

		hists = []
		max_h = 0
		for run in range(n_runs) :
			h, edges = np.histogram(pairs_overlap[run], bins = bin_edges)
			ax.plot(bin_centres, h)
			max_h = max(max_h, max(h))

		plt.xlabel('Number of overlapping gene mutations', fontsize = 24)
		plt.ylabel(r'frequency', fontsize = 28)
		plt.text(max_p, max_h, str(n_runs) + ' runs\n' + 'chromosome ' + chrom, 
			fontsize = 30, verticalalignment='top', horizontalalignment='right')

		outputFile = 'plots/' + 'pair_distribution_' + chrom + '_' + impact + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()

	def gene_pairs(self) :
		gene_pairs = []
		for item in file('gene_pairs_' + chrom + '_' + impact + '.txt') :
			item = item.split()
			gene_pairs.append(map(float, item))
	
		gene_list = []
		for item in file('genes_' + chrom + '_' + impact + '.txt') :
			gene_list.append(item.split()[0])

		print gene_pairs
		print gene_list

		gene_pairs = np.transpose(gene_pairs)

		fig = plt.figure(frameon=False, figsize=(10, 9), )
		fig.subplots_adjust(wspace = .0,top = .89, bottom = .2, left = .12, right = .95)
		ax = fig.add_subplot(111)

		cbar = plt.imshow(gene_pairs, interpolation='nearest', cmap = cm.nipy_spectral_r)
		cb = plt.colorbar(cbar)
		cb.set_label('N people with mutations on both genes')


		plt.xticks(range(len(gene_list)), gene_list, rotation='vertical', fontsize = 14,
			family = 'sans-serif')
		plt.yticks(range(len(gene_list)), gene_list, rotation='horizontal', fontsize = 14,
			family = 'sans-serif')

		text1 = str(len(gene_list)) + ' genes\n' + 'impact ' + \
			impact + '\n chromosome ' + chrom
		plt.text(len(gene_list)*.95, 0, text1, fontsize = 24, 
			verticalalignment='top', horizontalalignment='right')

		outputFile = 'plots/' + 'gene_pairs_' + chrom + '_' + impact + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()



if __name__ == '__main__':
	pt = PlotData()

	
	pt.patient_overlap()
	pt.gene_pairs()
