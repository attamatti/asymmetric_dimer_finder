#!/usr/bin/env python

#### update the path to chimera here ############

chimerapath='/fbs/emsoftware2/LINUX/fbsmi/Chimera-1.11.2-linux/bin/chimera'

#################################################

import sys
import subprocess
import os
import glob

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
    return(return_list)

def parse_loop(chunk):
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
        return(labels,data,False)
    else:
        for i in chunk[0]:
            if '_entity.id' in i:
                id = i.replace('/n','').split()[-1]
            if '_entity.pdbx_description' in i:
                i = i.replace('"',"'")
                name = i.replace('/n','').split("'")[1]
        labels = {id:0}
        data = [name]
        return({id:name},'X',True)

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
                line= i.split()
                i= i.replace('"',"'")
                splitline = i.split("'")[0].split()+[i.split("'")[1]]+i.split("'")[2].split()
            else:
                splitline = i.split()
               
            if splitline[labels['_entity.type']] == 'polymer':
                id = splitline[labels['_entity.id']]
                name = splitline[labels['_entity.pdbx_description']]
                entities[id] = name
        return(entities)
    elif skip == True:
        return(labels)
### make necessary dirs
if os.path.isdir('tmp') == False:
    subprocess.call(['mkdir','tmp'])
else:
    subprocess.call(['rm tmp/*'],shell=True)
tmpdir ='{0}/tmp'.format(os.getcwd())

if os.path.isdir('vis') == False:
    subprocess.call(['mkdir','vis'])


### read the mmcif file
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

if NMR_ensemble == False:    
    print('\nfile: {0}\tfound {1} polymer(s)'.format(sys.argv[1].split('.')[0],len(polyids)))
else:
    print('\nfile: {0}\tNMR emsemble\tfound {1} polymer(s)'.format(sys.argv[1].split('.')[0],len(polyids)))

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
        if chaindic[i][0] == chaindic[j][0] or chaindic[i][1] == chaindic[j][1]:
            if [j,i] not in matches and j != i:
                matches.append([i,j])
if len(matches) == 0:
    sys.exit('{0} matching chains'.format(len(matches)))

### write the chimera scripts to check for contacts
chimerascript = open('tmp/chimera_script.cmd','w')
chimerascript.write('open {0}/{1}\n'.format(os.getcwd(),sys.argv[1]))
for i in matches:
    if NMR_ensemble == False:
        chimerascript.write('findclash :.{0}&protein test  :.{1}&protein selectClashes true overlapCutoff -2.0 saveFile {3}/{2}_{0}{1}_A.txt\n'.format(i[0],i[1],sys.argv[1].split('.')[0],tmpdir))
        chimerascript.write('findclash :.{1}&protein test  :.{0}&protein selectClashes true overlapCutoff -2.0 saveFile {3}/{2}_{0}{1}_B.txt\n'.format(i[0],i[1],sys.argv[1].split('.')[0],tmpdir))
    elif NMR_ensemble == True:
        chimerascript.write('findclash #0.1:.{0}&protein test  #0.1:.{1}&protein selectClashes true overlapCutoff -2.0 saveFile {3}/{2}_{0}{1}_A.txt\n'.format(i[0],i[1],sys.argv[1].split('.')[0],tmpdir))
        chimerascript.write('findclash #0.1:.{1}&protein test  #0.1:.{0}&protein selectClashes true overlapCutoff -2.0 saveFile {3}/{2}_{0}{1}_B.txt\n'.format(i[0],i[1],sys.argv[1].split('.')[0],tmpdir))

chimerascript.close()

###run the chimera scripts
runchimera = subprocess.Popen('{0} --nogui {1}/chimera_script.cmd'.format(chimerapath,tmpdir), shell=True, stdout=subprocess.PIPE)
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
            #print('miss',[switched,switched2],i)
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
            #print('miss',[switched,switched2],i)
            misslistB.append([i,[switched,switched2]])
    print('{0} vs {1} {2}/{3} interactions found'.format(chain2, chain1,hitcount,len(d1))) 
    total = len(d1)+len(d2)
    if len(d1)+len(d2) > 0:
        return(float(len(hitlistA)+len(hitlistB))/float((len(d1)+len(d2))),misslistA,misslistB,hitlistA,hitlistB)


### read the chimera outputs
chifiles = glob.glob('{0}/*.txt'.format(tmpdir))
done = []

# determine which chains are forming higher level homo-oligomers
mmerdic = {}                    #{chain:# of other matches with > 0 contacts}
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
            mmerdic[filepair[0]] +=1
            mmerdic[filepair[1]] +=1

# screen output about oligomerization state
print('\nfinding homodimers\nchain\toligomer')
for i in chainkeys:
    if mmerdic[i] > 1:
        mmer = 'higher order oligomer (score: {})'.format(mmerdic[i])
    elif mmerdic[i] == 0:
        mmer = 'no oligomers detected'
    elif mmerdic[i] == 1:
        mmer = '** homodimer detected **'
    print('{0}\t{1}'.format(i,mmer))

mainout = open('asym_dimer-find.txt','a')
for i in chifiles:
    pair = (i.split('/')[-1].split('_')[1])
    if pair not in done and mmerdic[pair[0]] == 1:
        filename = '_'.join(i.split('_')[0:-1])
        afile = '{0}_A.txt'.format(filename)
        bfile = '{0}_B.txt'.format(filename)
        done.append(pair)
        contactdataa = open(afile,'r').readlines()
        contactdatab = open(bfile,'r').readlines()
        contactsa = parse_contact_data(contactdataa)
        contactsb = parse_contact_data(contactdatab)
        if len(contactsa) > 0 and len(contactsb) > 0:
            corr,missA,missB,hitA,hitB = match_contacts(contactsa,contactsb,pair[0],pair[1])
            output = open('{0}_colors.cmd'.format(filename.replace('/tmp/','/vis/'),pair[0],pair[1]),'w')
            print('{0} homodimer contact correlation = {1}\n'.format(pair,corr))
            
            output.write('~distance; ')
            output.write('transparency 0; ')
            output.write('color white; ')
            output.write('color cyan :.{0}; '.format(pair[0]))
            output.write('color hotpink :.{0}; '.format(pair[1]))
            for i in hitA:
                output.write('color green {0}; '.format(i))
            for i in missA:
                output.write('color red {0}; '.format(i[0].split()[0]))
                output.write('color blue {0}; '.format(i[0].split()[1]))
                output.write('color yellow {0}; '.format(i[1][0]))
                output.write('color orange {0}; '.format(i[1][1]))
            
            output.write('transparency 90 @/color=white\n')
            output.close()
            mainout.write('{0},{1},{2},{3}\n'.format(sys.argv[1],pair[0],pair[1],corr))
mainout.close()
