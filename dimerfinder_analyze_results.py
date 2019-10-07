#!/usr/bin/env python

import sys

errmsg = 'USAGE: dimerfinder_analyze_results.py <results file> <cc threshold> <rmsd threshold>'

try:
	data = open(sys.argv[1],'r').readlines()
except:
	sys.exit('ERROR reading results file\n{0}'.format(errmsg))
try:
	ccthresh = float(sys.argv[2])
except:
	sys.exit('ERROR cc thresh not specified\n{0}'.format(errmsg))
try:
	rmsdthresh = float(sys.argv[3])	
except:
	sys.exit('ERROR rmsd thresh not specified\n{0}'.format(errmsg))


allpdbs = {} 		#{file:filename, Chain1, Chain2, number of 1->2 contacts, number of 2->1 contacts, contact correlation,rmsd,min oligomer in file,protein name}

chain1,chain2,con12,con21,cc,rmsd,oligo,protein = range(0,8)

for line in data:
	i = line.split('%%')
	allpdbs[i[0]] = i[1:]

## find dimers 
count = 0
print('\npdb file\t\tcc\trmsd\tprotein')
for i in allpdbs:
	if allpdbs[i][oligo] == '2' and float(allpdbs[i][cc]) < ccthresh and float(allpdbs[i][rmsd]) < rmsdthresh:
		print('{0}#t\t{1:.2f}\t{2:.2f}\t{3}'.format(i,float(allpdbs[i][cc]),float(allpdbs[i][rmsd]),allpdbs[i][protein]))
		count+=1
if count == 0:
	print('** no matches found **')
