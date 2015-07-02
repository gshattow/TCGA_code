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
datafile = '50000_variant_effect_output.txt'


class ReadData :
	def read_vep_file(self, datafile) :
		genes_all = []
		patients_all = []

		location_index = 0
		extra_index = 0
		gene_index = 0
		patient_index = 0
		genes_all = []
		patients_all = []
		chromosomes = []
		
		ii = 0
		tic = time.clock()
		for item in file(datafile) :
			ii = ii + 1
			if item.startswith('##'):
#				print ii
				continue
			elif item.startswith('#') :
				print 'line number', ii, 'has the index names'
				item = item.split()
				location_index = item.index('Location')
				extra_index = item.index('Extra')
#			if not item.startswith('#') :
			else :
#				print ii
				item = item.split()
				loc = item[location_index].split(':')
#				chromosomes.append(loc[0])
				extra = item[extra_index].replace(';',' ').replace('=', ' ').split()
#				print extra
				if set(['SYMBOL', 'IND']) <= set(extra) :
#				if 'SYMBOL' and 'IND' in extra :
# 					print ii, len(extra), extra
					if extra[extra.index('SYMBOL_SOURCE') + 1] == 'HGNC' :
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
		
		print len(genes_all)
		print len(patients_all)

		n_g = len(genes_all)
		n_p = len(patients_all)
		
		print genes_all
		print patients_all

############################################################



if __name__ == '__main__':
	rd = ReadData()
	rd.read_vep_file(datafile)
