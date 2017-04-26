from __future__ import print_function, division

#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Copyright (c) 2016 Wei-Xin Liu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import os, sys, gzip
import pandas as pd
from itertools import dropwhile, islice
from pandas import DataFrame, Series
from time import sleep

def gffClient(filename, selected_fld=None, attr_fld=None, LEN=None, To_file=False):
	'''
	Parameters:
		filename: the gff file, either ends with `gff.gz` or `gff`
		selected_fld: int LIST. the fileds indexes (0-based, 0~8) that wanted to extract. There 9 fileds + 1 [comments]:
			<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
		attr_fld: STRING LIST. the fileds in [attributes] that want to extract. It can be the <gene_id>, <transcript_id>, <gene_name>, ...
		LEN: the number of items in the file named `filename`
		To_file: whether save the result to file derectly
		
	Returns: 
		Data: int 0 or DataFrame. 
	'''
	header = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
	if filename.endswith('gff.gz'):
		gffopen = gzip.open
		outfile = '.'.join(filename.split('.')[:-2])
	elif filename.endswith('gff'):
		gffopen = open
		outfile = '.'.join(filename.split('.')[:-1])
	else:
		raise IOError, r'the format of the file are not gff!'
		
	with gffopen(filename, 'rb') as handle:
		idx = 0
		if LEN < 10000:
			LEN = None
		if To_file:
			OUT = gzip.open('{}.tsv.gz'.format(outfile), 'wb')
			observe = '\t'.join([header[selected_idx] for selected_idx in selected_fld] + attr_fld)
			OUT.write(observe)
			OUT.write('\n')
			
			#for line in dropwhile(lambda line: line.startswith('#'), handle):
			for line in islice(handle, 0, None):
				if LEN and (idx % 10000 == 0):
					progbar(LEN, idx)
				
				line = line.strip().split('\t')
				try:
					attr = [itm.strip().split('=') for itm in line[8].split(';')]
				except Exception, e:
					OUT.write(line[0])
					OUT.write('\n')
					continue
				attr = {itm[0].strip():itm[1] for itm in attr if len(itm) > 1}
				
				attr_str = []
				for j in attr_fld:
					try:
						attr_str.append(attr[j] )
					except Exception, e:
						#print(e)
						attr_str.append('.')

				observe = '\t'.join([line[selected_idx] for selected_idx in selected_fld] + attr_str)
				OUT.write(observe)
				OUT.write('\n')
				
				idx = idx + 1
				
			#print('\n')	
			OUT.close()	
			return 0
				
		else:
			Data = DataFrame()
			for line in dropwhile(lambda line: line.startswith('#'), handle):
				if LEN and (idx % 10000 == 0):
					progbar(LEN, idx)
				
				line = line.strip().split('\t')
				attr = [itm.strip().split('=') for itm in line[-1].split(';')]
				attr = {itm[0].strip():itm[1] for itm in attr if len(itm) > 1}
			
				fileds = {}
				for selected_idx in selected_fld:
					fileds[header[selected_idx] ] = line[selected_idx]
				for j in attr_fld:
					try:
						fileds[j] = attr[j]
					except Exception, e:
						#print(e)
						fileds[j] = '.'
				
				oneItem = Series(fileds)
				if idx == 0:
					oneItem.name = str(idx)
					Data = DataFrame(oneItem)
				else:
					Data[str(idx)] = oneItem
							
				idx = idx + 1
			
			print('\n')	
			Data = Data.T
			return Data
	
def progbar(LEN, i):
	Interval = LEN // 60
	if (i) % Interval == 0:
		j = '#' * ((i+1) // Interval)
		sys.stdout.write(str(int(((i+1)/LEN)*100))+'% ||'+j+'->'+'\r')
		sys.stdout.flush()
		
def work():
	# ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
	DIR = r'/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_Mus/2+.Count2FPKM/GRCm38.p5'
	data = gffClient(r'/'.join([DIR, 'GCF_000001635.25_GRCm38.p5_genomic.gff.gz'] ), selected_fld=[0, 1, 2, 3, 4, 6], attr_fld=['Dbxref', 'gene', 'transcript_id', 'gene_synonym'], LEN=2286189, To_file=True)
	if isinstance(data, DataFrame):
		data['length'] = pd.to_numeric(data['end']) - pd.to_numeric(data['start'])
		del data['end']
		del data['start']
		data.to_csv(r'/'.join([DIR, 'GCF_000001635.25_GRCm38.p5_genomic.csv'] ) )
			
if __name__ == '__main__':
	#work()
	pass

	
