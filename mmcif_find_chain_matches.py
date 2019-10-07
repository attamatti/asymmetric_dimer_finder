#!/usr/bin/env python

#### update the path to chimera here ############

chimerapath='/fbs/emsoftware2/LINUX/fbsmi/Chimera-1.11.2-linux/bin/chimera'

#################################################

import sys
import subprocess
import os
import glob
FNULL = open(os.devnull, 'w')

def parse_mmcif(ciffile):
    cifdata = open(ciffile,'r').readlines()
    chunks = []
    n=-1
    indata= False
    for i in cifdata:
        if '#' in i:
            indata=True
            n+=1    
            chunks.append([])
        if indata==True:
            chunks[n].append(i)
    return(chunks)

def return_chunks(data,search_string):
    return_list = []
    for i in data:
        for j in i:
            if search_string in j:
                return_list.append(i)
    if len(return_list) > 0:
        return(return_list)
    else:
        stats.write('{0},{1},{2}\n'.format(sys.argv[1].split('/')[-1],'NOHEADER',search_string))
        sys.exit('FILE ERROR: no {0} header found in mmcif file'.format(search_string))

def parse_loop(chunk):
	try:
		if 'loop_\n' in chunk[0]: 
			labels = {}         #{label:column_number}
			indata = False
			data = []
			trip = False
			n=0
			inww = False
			for i in chunk[0][2:]:
				if i[0] != '_':
					indata= True 
				if indata == False:
					labels[i.replace('\n','').strip()] = n
				elif indata == True:   
					if inww == True:
						data[-1] = data[-1]+[i]
						if i[0] == ';':
							inww = not inww
							trip = True
					if inww ==False and trip != True:
						if i[0] ==';':
							data[-1] = data[-1]+[i]
							inww = not inww
						elif chunk[0][n+1][0]== ';':
							data[-1] = data[-1]+[i.replace('\n','')]
						else:
							data.append([i.replace('\n','')])
					trip=False
				n+=1
			n=0
			for i in data:
				if len(i) > 1:
					data[n] = ['  '.join(i).replace('\n','').replace(';','')]
					n+=1
			return(labels,data,False)
		else:
			for i in chunk[0]:
				if '_entity.id' in i:
					id = i.replace('/n','').split()[-1]
				if '_entity.pdbx_description' in i:
					if '"' in i or "'" in i:
						i = i.replace('"',"'")
						name = i.replace('/n','').split("'")[1]
					else:
						name = i.split()[1]
			labels = {id:0}
			data = [name]
			return({id:name},'X',True)
	except:
		stats.write('{0},FILE_READ_ERROR,parse_loop'.format(sys.argv[1].split('/')[-1]))
		sys.exit('FILE READ ERROR in parse_loop - skipping file')
def parse_loop_strand(chunk):
    labels = {}         #{label:column_number}
    indata = False
    data = []
    n=0
    inww = False
    comblines = []
    for i in chunk[0][2:]:
        if i[0] != '_':
            indata= True 
        
        if indata == False:
            labels[i.replace('\n','').strip()] = n
        
        elif indata == True:
            comblines.append(i)

        n+=1
    splitdat = []
    comblines = [''.join(comblines)]
    comblines[0] = comblines[0].replace('\n',' ') 
    for j in comblines[0].split():
        if inww == False:
            splitdat.append(j)
        if inww == True:
            splitdat[-1] = splitdat[-1]+j
        if j[0] == ';':
            inww = not inww
    findat = [splitdat[i:i + len(labels)] for i in xrange(0, len(splitdat), len(labels))]
    return(labels,findat)

def parse_single(chunk):
    chainsdic = {}           #{chain: EntityID}
    for i in chunk[0]:
        if i.split()[0] == '_entity_poly.pdbx_strand_id':
            chains = i.split()[1]
        if i.split()[0] == '_entity_poly.entity_id':
            eid = i.split()[1]
    for i in chains.split(','):
        chainsdic[i] = eid
    return(chainsdic)

def get_polys(labels,data,skip):
    if skip == False:
        polys = []
        entities = {}           #{entityID:name}
        for i in data:
            i = i[0]
            if "'" in i or '"' in i:
                ccount1 = 0
                for j in i.split():
                    if "'" in j:
                        ccount1+=1 
                if ccount1 ==1:
                    splitline = i.split()
                else:
                    line= i.split()
                    i= i.replace("'",'"')
                    splitline = i.split('"')[0].split()+[i.split('"')[1]]+i.split('"')[2].split()
            else:
                splitline = i.split()
            try:
                if splitline[labels['_entity.type']] == 'polymer':
                    id = splitline[labels['_entity.id']]
                    name = splitline[labels['_entity.pdbx_description']]
                    entities[id] = name
            except:
                print('bad input line in get_polys: skipping line {0}'.format(splitline))
        return(entities)
    elif skip == True:
        return(labels)

def parse_assembly(chunk,ss):	
	aidlist = []	
	oligolist = []
	dimer = False
	if 'loop_\n' in chunk[0]:
		aidlist = []		
		labels = {}         #{label:column_number}
		indata = False
		data = []
		n=0
		inww = False	
		comblines = []
		for i in chunk[0][2:]:
			if i[0] != '_':
				indata= True 
	        
			if indata == False:
				labels[i.replace('\n','').strip()] = n
	        
			elif indata == True:
				comblines.append(i)
	
			n+=1
		splitdat = []
		comblines = [''.join(comblines)]
		comblines[0] = comblines[0].replace('\n',' ')		
		
		for j in comblines[0].split():
			if inww == False:
				splitdat.append(j)
			if inww == True:
				splitdat[-1] = splitdat[-1]+j
			if j[0] == "'" or j[-1] == "'":
				inww = not inww
		findat = [splitdat[i:i + len(labels)] for i in xrange(0, len(splitdat), len(labels))]
		for i in findat:
			oligolist.append(i[labels[ss]])
			try:			
				aidlist.append(i[labels['_pdbx_struct_assembly.id']])
			except:
				aidlist.append(i[labels['_pdbx_struct_assembly_gen.assembly_id']])
	else:
		for i in chunk[0]:		
			if ss in i:
				oligolist = [i.replace('/n','').split()[-1]]
		aidlist = [1]	
	return(oligolist,aidlist)


### make necessary dirs
if os.path.isdir('tmp') == False:
    subprocess.call(['mkdir','tmp'])

subprocess.call(['rm tmp/*'],shell=True,stderr=FNULL)
tmpdir ='{0}/tmp'.format(os.getcwd())

if os.path.isdir('vis') == False:
    subprocess.call(['mkdir','vis'])
subprocess.call(['touch','stats.txt'])
stats =open('stats.txt','a')

### read the mmcif file
print('\ninput path: {0}'.format(sys.argv[1]))
pdbid = sys.argv[1].split('/')[-1].split('.')[0]
chunks = parse_mmcif(sys.argv[1])

# determine if it's a NMR ensemble - needs to be treated differently
NMR_ensemble = False
for i in chunks:
    for j in i:
        if 'pdbx_nmr' in j:
            NMR_ensemble = True
            break

entity_chunk = (return_chunks(chunks,'_entity.id'))
elabels,edata,skip = parse_loop(entity_chunk)
polyids = (get_polys(elabels,edata,skip))            #{entityID:name}
strand_chunk = (return_chunks(chunks,'_entity_poly.pdbx_strand_id'))


#try:
dimer_assembly_list,dimerids = parse_assembly(return_chunks(chunks,'_pdbx_struct_assembly.oligomeric_details'),'_pdbx_struct_assembly.oligomeric_count')
oligo_assembly_list,oligoids = parse_assembly(return_chunks(chunks,'_pdbx_struct_assembly_gen.asym_id_list'),'_pdbx_struct_assembly_gen.asym_id_list')
#except:
#	dimer_assembly_list = [0]	
#	oligo_assembly_list = []
#print(dimer_assembly_list)
#print(dimerids)
oligodic = dict(zip(dimerids,dimer_assembly_list))		# {assembly ID:oligomeric state }
#print ('odic',oligodic)
#print(oligo_assembly_list)
#print(oligoids)
assemblydic = dict(zip(oligo_assembly_list,oligoids))	# {assembly chains:AssemblyID}
#print('adic',assemblydic)

## get the entiy IDs and protein names of the different chains
chaindic = {}           #{chain:entityID}
if strand_chunk[0][1] == 'loop_\n':
    slabels,sdata = parse_loop_strand(strand_chunk)
    for i in sdata:
        if ',' in i[slabels['_entity_poly.pdbx_strand_id']]:
            for chain in i[slabels['_entity_poly.pdbx_strand_id']].split(','):
                chaindic[chain] = i[slabels['_entity_poly.entity_id']]
        else:
            chaindic[i[slabels['_entity_poly.pdbx_strand_id']]] = i[slabels['_entity_poly.entity_id']]
else:
    chaindic = parse_single(strand_chunk)

## treat NMR ensembles differently because ofthe way the models are structured
## use model #0.1/#1.1 instead of #0/#1
if NMR_ensemble == False:    
    print('file: {0}\tfound {1} polymer(s)'.format(sys.argv[1].split('.')[0],len(polyids)))
else:
    print('file: {0}\tNMR emsemble\tfound {1} polymer(s)'.format(sys.argv[1].split('.')[0],len(polyids)))
if len(polyids) == 0:
    stats.write('{0},{1}\n'.format(sys.argv[1].split('/')[-1],'NO_PROTEIN_POLYS'))
    sys.exit('FILE ERROR: no polymers found')

print('Chain\tID\tName')


### identify matching chains
chainkeys = list(chaindic)
chainkeys.sort()
for i in chainkeys:                      #coverted to {chain:[entity_id,Name]}
    print('{0}\t{1}\t{2}'.format(i,chaindic[i],polyids[chaindic[i]]))
    chaindic[i] = [chaindic[i],polyids[chaindic[i]]]
matches = []
for i in chaindic:
    for j in chaindic:
        if chaindic[i][0] == chaindic[j][0]:
            if [j,i] not in matches and j != i:
                matches.append([i,j])
if len(matches) == 0:
    stats.write('{0},{1}\n'.format(sys.argv[1].split('/')[-1],'NO_IDENTICAL_CHAINS'))
    sys.exit('{0} matching chains'.format(len(matches)))
if len(matches) > 20:
	stats.write('{0},TOO_MANY,{1}\n'.format(sys.argv[1].split('/')[-1],len(matches)))	
	sys.exit('Found {0} matching chains; too many matching chains - skipping'.format(len(matches)))

### write the chimera scripts to check for contacts
chimerascript = open('tmp/chimera_script.cmd','w')
chimerascript.write('open {0}/{1}\n'.format(os.getcwd(),sys.argv[1]))
for i in matches:
    if NMR_ensemble == False:
        chimerascript.write('findclash :.{0}&protein test  :.{1}&protein selectClashes true overlapCutoff -1.0 saveFile {3}/{2}_{0}{1}_A.txt\n'.format(i[0],i[1],pdbid,tmpdir))
        chimerascript.write('findclash :.{1}&protein test  :.{0}&protein selectClashes true overlapCutoff -1.0 saveFile {3}/{2}_{0}{1}_B.txt\n'.format(i[0],i[1],pdbid,tmpdir))
    elif NMR_ensemble == True:
        chimerascript.write('findclash #0.1:.{0}&protein test  #0.1:.{1}&protein selectClashes true overlapCutoff -1.0 saveFile {3}/{2}_{0}{1}_A.txt\n'.format(i[0],i[1],pdbid,tmpdir))
        chimerascript.write('findclash #0.1:.{1}&protein test  #0.1:.{0}&protein selectClashes true overlapCutoff -1.0 saveFile {3}/{2}_{0}{1}_B.txt\n'.format(i[0],i[1],pdbid,tmpdir))

chimerascript.close()

###run the chimera scripts
runchimera = subprocess.Popen('{0} --nogui {1}/chimera_script.cmd'.format(chimerapath,tmpdir), shell=True, stdout=subprocess.PIPE,stderr=FNULL)
chimeraout = runchimera.stdout.read()


def parse_contact_data(cdata):
    contactslist = []        
    for i in cdata:
        i = i.replace('#0.1:',':')
        if i[0] == ':':
            c1 = i.split()[0]
            c2 = i.split()[1]
            contactslist.append('{0} {1}'.format(c1,c2))
    return(contactslist)

def match_contacts(d1,d2,chain1,chain2):
    hitcount = 0
    d11 = [x.split()[0] for x in d1]
    d22 = [x.split()[0] for x in d2]
    misslistA,misslistB = [],[]
    hitlistA,hitlistB = [],[]
    for i in d1:
        switched = '.{0}@'.format(chain2).join([i.split()[0].split('.')[0],i.split()[0].split('@')[-1]])
        switched2 = '.{0}@'.format(chain1).join([i.split()[1].split('.')[0],i.split()[1].split('@')[-1]])
        if switched in d22:
            hitcount+=1
            hitlistA.append(i)
        else:
            misslistA.append([i,[switched,switched2]])
    print('\n{0} vs {1} {2}/{3} interactions found'.format(chain1,chain2,hitcount,len(d1))) 
    hitcount = 0
    for i in d2:
        switched = '.{0}@'.format(chain1).join([i.split()[0].split('.')[0],i.split()[0].split('@')[-1]])
        switched2 = '.{0}@'.format(chain2).join([i.split()[1].split('.')[0],i.split()[1].split('@')[-1]])
        if switched in d11:
            hitcount+=1
            hitlistB.append(i)
        else:
            misslistB.append([i,[switched,switched2]])
    print('{0} vs {1} {2}/{3} interactions found'.format(chain2, chain1,hitcount,len(d1))) 
    total = len(d1)+len(d2)
    if len(d1)+len(d2) > 0:
        return(float(len(hitlistA)+len(hitlistB))/float((len(d1)+len(d2))),misslistA,misslistB,hitlistA,hitlistB)


### read the chimera outputs
chifiles = glob.glob('{0}/{1}*.txt'.format(tmpdir,pdbid))
done = []

# determine which chains are forming higher level homo-oligomers
mmerdic = {}                    #{chain:# of other matches with > 0 contacts}
partners = {}					#{chain:[partner1,partner2]}
for i in chifiles:
	filehalf = i.split('_')[-1].replace('.txt','')
	filepair = [i.split('_')[-2][0],i.split('_')[-2][1]]
	try:
		float(mmerdic[filepair[0]])
	except:
		mmerdic[filepair[0]] = 0
	try:
		float(mmerdic[filepair[1]])
	except:
		mmerdic[filepair[1]] = 0
	if filehalf == 'A':
		contacts = parse_contact_data(open(i,'r').readlines())
		if len(contacts) > 0:
			try:
				partners[filepair[0]].append(filepair[1])
			except:
				partners[filepair[0]] = [filepair[1]]
			try:
				partners[filepair[1]].append(filepair[0])
			except:
				partners[filepair[1]] = [filepair[0]]
			mmerdic[filepair[0]] +=1
			mmerdic[filepair[1]] +=1
		elif len(contacts) == 0:
			partners[filepair[0]] = []
			partners[filepair[1]] = []


## you've hooked up with everyone your partner has hooked up with :P
partnerscores = {}				#{chain:total score}
for i in partners:
	partnerscores[i] = sum([mmerdic[x] for x in partners[i]])

# screen output about oligomerization state
print('\nfinding homodimers\nchain\toligomer')

for i in chainkeys:
	try:
		if partnerscores[i] > 1:
			mmer = 'higher order oligomer (score: {0})'.format(partnerscores[i])
		elif partnerscores[i] == 0:
			mmer = 'no oligomers detected'
		elif partnerscores[i] == 1:
			if partnerscores[partners[i][0]] == 1:			
				mmer = '** homodimer detected with chain {0}**'.format(partners[i][0])
			else:
				mmer = 'homodimer detected with chain {0} but not reciprocated'.format(partners[i][0]) 		
		print('{0}\t{1}'.format(i,mmer))
	except:
		print('{0}\tskipping chain - not protein'.format(i))

## compare the results of the chimera scripts and write main output
mainout = open('asym_dimer-find.txt','a')
homodimers = False
hdcount = 0
reds,blues,yellows,oranges = [],[],[],[]
for i in chifiles:
	homodimer = False
	pair = (i.split('/')[-1].split('_')[1])
	if pair not in done and partnerscores[pair[0]] == 1 and partnerscores[pair[1]] == 1:
		filename = '_'.join(i.split('/')[-1].split('_')[0:-1])
		afile = 'tmp/{0}_A.txt'.format(filename)
		bfile = 'tmp/{0}_B.txt'.format(filename)
		done.append(pair)
		contactdataa = open(afile,'r').readlines()
		contactdatab = open(bfile,'r').readlines()
		contactsa = parse_contact_data(contactdataa)
		contactsb = parse_contact_data(contactdatab)

		if len(contactsa) > 0 and len(contactsb) > 0:
			homodimer = True
			corr,missA,missB,hitA,hitB = match_contacts(contactsa,contactsb,pair[0],pair[1])
			print('{0} homodimer contact correlation = {1}'.format(pair,corr))
            
	    	### calculate the RMSD of the chains that were
			chimerarmsd = open('{0}/chirmsd.cmd'.format(tmpdir),'w')
			if NMR_ensemble ==True:
				chimerarmsd.write('open #0 {3}/{0};open #1 {3}/{0};mmaker #1.1:.{1} #0.1:.{2}'.format(sys.argv[1],pair[0],pair[1],os.getcwd()))
			else:
				chimerarmsd.write('open #0 {3}/{0};open #1 {3}/{0};mmaker #1:.{1} #0:.{2}'.format(sys.argv[1],pair[0],pair[1],os.getcwd()))
			chimerarmsd.close()
	    
	    	### run script 
 			runchimera2 = subprocess.Popen('{0} --nogui {1}/chirmsd.cmd'.format(chimerapath,tmpdir), shell=True, stdout=subprocess.PIPE,stderr=FNULL)
			chimeraout2 = runchimera2.stdout.read()	    
	    		
			for i in chimeraout2.split('\n'):
				if 'RMSD' in i:
					rmsd = float(i.split()[-1].strip(')'))
					print('RMSD between chains {0}/{1} = {2}'.format(pair[0],pair[1],rmsd))			
			showlist = []			
			output = open('vis/{0}_colors.cmd'.format(filename.replace('/tmp/','/vis/'),pair[0],pair[1]),'w')
			output.write('transparency 0; ')
			output.write('color white; ')
			output.write('color cyan :.{0}; '.format(pair[0]))
			output.write('color hotpink :.{0}; '.format(pair[1]))
			hitcolor = []			
			for i in hitA:
				hitcolor.append('{0}'.format(i))
				showlist.append(i.split('@')[0])
				showlist.append(i.split()[-1].split('@')[0])
			for i in missA:
				reds.append(i[0].split()[0])
				showlist.append(i[0].split()[0].split('@')[0])
				blues.append(i[0].split()[1])
				showlist.append(i[0].split()[1].split('@')[0])				
				yellows.append(i[1][0])
				showlist.append(i[1][0].split('@')[0])				
				oranges.append(i[1][1])	
				showlist.append(i[1][1].split('@')[0])		
			output.write('show {0}\n'.format(' '.join(showlist)))			
			output.write('transparency 90 @/color=white\n')
			if len(hitA) > 0:			
				output.write('color green {0}; '.format(' '.join(hitcolor)))			
			if len(reds) > 0:
				output.write('color red {0}; '.format(' '.join(reds)))
			if len(blues) > 0:
				output.write('color blue {0}; '.format(' '.join(blues)))
			if len(yellows) > 0:
				output.write('color yellow {0}; '.format(' '.join(yellows)))
			if len(oranges) > 0:
				output.write('color orange {0}; '.format(' '.join(oranges)))
			output.close()
		
		### finally check to see if dimer is designated as biological assembly in mmcif file
			bdlist = []
			for i in assemblydic:
				if pair[0] in i and pair[1] in i:
					bdlist.append(assemblydic[i])
			hitlist = []
			for i in bdlist:
				hitlist.append(oligodic[i])
			for i in bdlist:
				if oligodic[i] == min(hitlist):
					biodimer = oligodic[i]
			if bdlist == []:
				biodimer = 1
			print('smallest designated biological oligomer: {0}-mer'.format(biodimer))
	



		## write the main output deliniated with %% to deal with spaces,quotes,naostandard characters and other bullshit
		## format is filename, Chain1, Chain2, number of 1->2 contacts, number of 2->1 contacts, contact correlation,rmsd,min biological dimer(T/F),protein name  
			mainout.write('{0}%%{1}%%{2}%%{3}%%{4}%%{5}%%{6}%%{7}%%{8}\n'.format(sys.argv[1],pair[0],pair[1],len(contactsa),len(contactsb),corr,rmsd,biodimer,chaindic[pair[0]][1]))
			if homodimer == True:
				homodimers +=1 
if homodimers >  0:
	stats.write('{0},{1},{2}\n'.format(sys.argv[1].split('/')[-1],'HOMODIMER(s)',homodimers))
else:
	stats.write('{0},{1}\n'.format(sys.argv[1].split('/')[-1],'NOHOMODIMERS'))

mainout.close()
