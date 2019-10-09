#!/usr/bin/env python

import sys
import os
import subprocess
FNULL = open(os.devnull, 'w')

#### update the path to chimera here ############

chimerapath='/fbs/emsoftware2/env-modules/install/chimera/1.13.1/bin/chimera'

#################################################
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


cwd = os.getcwd()
## find dimers 
hits = 0
hitlist = []
rotang = 0
print('\npdb file\tcc\trmsd\t#int\trot\tseqid\tprotein')

for i in allpdbs:
	print(i)
	if allpdbs[i][oligo] == '2' and float(allpdbs[i][cc]) < ccthresh and float(allpdbs[i][rmsd]) < rmsdthresh:
		hitlist.append(i)

## write chimera script to find the symmetry relations 
		chiscript = open('tmp/chirot.cmd','w')
		chiscript.write('open {0}/{1}; open {0}/{1};mmaker #1:.{2} #0:.{3} verbose True;measure rotation #0 #1'.format(cwd,i,allpdbs[i][chain1],allpdbs[i][chain2]))
		chiscript.close()		
		runchimera2 = subprocess.Popen('{0} --nogui {1}/chirot.cmd'.format(chimerapath,'{0}/tmp'.format(os.getcwd())), shell=True, stdout=subprocess.PIPE,stderr=FNULL)
		chimeraout2 = runchimera2.stdout.read()	    
		con = 0
		for j in chimeraout2.split('\n'):
			if 'Rotation angle (degrees)' in j:
				rotang = float(j.split()[-1])
			if 'Residues:' in j:
				seq1 = chimeraout2.split('\n')[con+2]
				seq2 = chimeraout2.split('\n')[con+4]
			con+=1
		allpdbs[i].append(rotang)
	#calculate sequence overlap the the supposedly 'identical' halves
		seqreplace = [x.replace('#0','#1').replace('.{0}'.format(allpdbs[i][chain2]),'.{0}'.format(allpdbs[i][chain1])) for x in seq2.split(',')]
		seqid = len(set(seq1.split(','))&set(seqreplace))/float(len(seq1.split(',')))
		allpdbs[i].append(seqid)
# now filter on more characteristics:
finalhits = []
for i in hitlist:
	if 10 < allpdbs[i][-2] < 170 and allpdbs[i][-1] > 0.90:
		finalhits.append([i,str(round(float(allpdbs[i][cc]),2)),str(round(float(allpdbs[i][rmsd]),2)),str(allpdbs[i][con21]),str(round(allpdbs[i][-2],2)),str(allpdbs[i][-1]),allpdbs[i][protein].replace('\n','')])


finalhits.sort(key=lambda x:float(x[3]),reverse=True)
if len(finalhits) > 0:
	for i in finalhits:
		print'\t'.join(i)
else:
	print('** no matches found **')
