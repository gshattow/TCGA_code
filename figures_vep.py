import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import params
import pylab
import matplotlib.pyplot as plt
from random import sample
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable


font = {'family':'serif',
	'size':14	}
plt.rc('font', **font)

plot_dir = 'plots/'
OutputFormat = '.png'


SIFT = "%03d" % int(params.sift*100)



class ReadData :
	def read_pairs_overlap(self) :
		file = 'patient_pairs_overlap_c' + params.chrom + '_s' + SIFT + '.txt'
		pairs_overlap = np.loadtxt(file, unpack=True)
		pairs_overlap = np.transpose(pairs_overlap)
		
		return pairs_overlap

	def read_gene_pairs(self) :
		file = 'gene_pairs_c' + params.chrom + '_s' + SIFT + '.txt'
		gene_pairs = np.loadtxt(file, unpack=True)

		return gene_pairs

	def read_genes(self) :
		file = 'genes_c' + params.chrom + '.txt'
		genes_list = np.loadtxt(file, dtype=np.str)
		
		return genes_list


class PlotData :		

	def patient_overlap(self, pairs_overlap) :
		
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
		text1 = str(len(gene_list)) + ' genes\n' + 'SIFT < ' + \
			str(params.sift) + '\n chromosome ' + params.chrom + '\n' + \
			str(params.n_runs) + ' runs'
		plt.text(.95, .95, text1, fontsize = 24, 
			verticalalignment='top', horizontalalignment='right',
			transform = ax.transAxes)

		outputFile = plot_dir + 'pair_distribution_c' + params.chrom + '_s' + \
			SIFT + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()

	def gene_pairs(self, gene_pairs) :


		fig = plt.figure(frameon=False, figsize=(10, 9), )
		fig.subplots_adjust(wspace = .0,top = .99, bottom = .15, left = .15, right = .9)
		ax = fig.add_subplot(111)

		plt.xticks(range(len(gene_list)), gene_list, rotation='vertical', fontsize = 10,
			family = 'sans-serif')
		plt.yticks(range(len(gene_list)), gene_list, rotation='horizontal', fontsize = 10,
			family = 'sans-serif')

		text1 = str(len(gene_list)) + ' genes\n' + 'SIFT < ' + \
			str(params.sift) + '\n chromosome ' + params.chrom 
		plt.text(0.95, 0.95, text1, fontsize = 24, 
			verticalalignment='top', horizontalalignment='right',
			transform = ax.transAxes)

		cbar = plt.imshow(gene_pairs, interpolation='nearest', cmap = cm.nipy_spectral_r)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="2%", pad=0.05)
		cb = plt.colorbar(cbar, cax = cax)
		cb.set_label('N people with mutations on both genes', fontsize = 14)

		outputFile = plot_dir + 'gene_pairs_c' + params.chrom + '_s' + SIFT + OutputFormat
		plt.savefig(outputFile)  
		print 'Saved file to', outputFile
		plt.close()



if __name__ == '__main__':
	rd = ReadData()
	pd = PlotData()

	pairs_overlap = rd.read_pairs_overlap()
	gene_pairs = rd.read_gene_pairs()
	gene_list = rd.read_genes()
	
	pd.patient_overlap(pairs_overlap)
	pd.gene_pairs(gene_pairs)
