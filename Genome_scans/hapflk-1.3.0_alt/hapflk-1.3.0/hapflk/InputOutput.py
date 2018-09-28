import sys
import data
import re
import numpy as np
from bisect import insort
from hapflk import missing
try:
    import argparse
    #### create a common parser for IO operations
    def plink_file(string):
        try:
            open(string+'.ped')
        except IOError:
            msg="can't open PED file %s"%string+'.ped'
            raise argparse.ArgumentTypeError(msg)
        try:
            open(string+'.map')
        except IOError:
            msg="can't open MAP file %s"%string+'.map'
            raise argparse.ArgumentTypeError(msg)
        return string
    def plink_bfile(string):
        try:
            open(string+'.bed')
        except IOError:
            msg="can't open BED file %s"%string+'.bed'
            raise argparse.ArgumentTypeError(msg)
        try:
            open(string+'.fam')
        except IOError:
            msg="can't open FAM file %s"%string+'.fam'
            raise argparse.ArgumentTypeError(msg)
        try:
            open(string+'.bim')
        except IOError:
            msg="can't open BIM file %s"%string+'.bim'
            raise argparse.ArgumentTypeError(msg)
        return string
    io_parser=argparse.ArgumentParser(add_help=False,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    input_group=io_parser.add_argument_group('Input Files','Available formats')
    input_group.add_argument('--ped',metavar='FILE',help='PED file')
    input_group.add_argument('--map',metavar='FILE',help='MAP file')
    input_group.add_argument('--file',metavar='PREFIX',help='PLINK file prefix (ped,map)',type=plink_file)
    input_group.add_argument('--bfile',metavar='PREFIX',help='PLINK bfile prefix (bim,fam,bed)',type=plink_bfile)
    param_group=io_parser.add_argument_group('Data format options','')
    param_group.add_argument('--miss_geno',help='Missing Genotype Code',metavar='C',default='0')
    param_group.add_argument('--miss_pheno',help='Missing Phenotype Code',metavar='C',default='-9')
    param_group.add_argument('--miss_parent',help='Missing Parent Code',metavar='C',default='0')
    param_group.add_argument('--miss_sex',help='Missing Sex Code',metavar='C',default='0')
    map_group=io_parser.add_argument_group('SNP selection','Filter SNP by genome position')
    map_group.add_argument('--chr',help='Select chromosome C',metavar='C')
    map_group.add_argument('--from',dest='posleft',help='Select SNPs with position > x (in bp/cM) Warning : does not work with BED files',metavar='x',default=0,type=float)
    map_group.add_argument('--to',dest='posright',help='Select SNPs with position < x (in bp/cM) Warning : does not work with BED files', metavar='x',default=np.inf,type=float)
    map_group.add_argument('--other_map',help='Use alternative map (genetic/RH)',action='store_true',default=False)
except:
    pass

defaultParams= {
    "MissGeno":'0',
    "MissPheno":'-9',
    "MissParent":'0',
    "MissSex":'0'
    }

key="kshuwewekiuvf"

def _get_snpIdx(myMap,chrom,left,right,othermap=False):
    '''
    Return index of SNPs in map

    Parameters
    --------------
    myMap : the map
    chrom : chromosome
    left : leftmost position
    right : rightmost position
    othermap : use other (genetic / RH) map scale
    '''
    if othermap:
        getMap=myMap.genetMap
    else:
        getMap=myMap.physMap
    mySnps=[ x for x in getMap(chrom,xleft=left,xright=right)]
    mySnpIdx=[]
    for i,s in enumerate(myMap.input_order):
        if s in mySnps:
            mySnpIdx.append(i)
    return mySnpIdx

##### generic input function ####
def parseInput(options,params=defaultParams):
    '''
    Read input files 
    
    Parameters
    --------------
    options : options from io_parser
    params : parameters for coding used in input files, overidden using options

    Returns
    ----------
    res: dictionary with keys:
         -- dataset : an instance from the hapflk Dataset class
         -- map : an instance from the hapflk Map class

    See also
    ----------
    data.Dataset
    data.Map
    '''

    ## get params from options
    if options.miss_geno:
        params["MissGeno"]=options.miss_geno
    if options.miss_pheno:
        params["MissPheno"]=options.miss_pheno
    if options.miss_parent:
        params["MissParent"]=options.miss_parent
    if options.miss_sex:
        params["MissSex"]=options.miss_sex
    ## read input
    ## 1. case with separate ped / map : ped needs a header line
    if options.ped:
       if options.map:
           ### get map from file and list of SNPs to keep
           myMap=parseMapFile(options.map)
           if options.other_map:
               getMap=myMap.genetMap
           else:
               getMap=myMap.physMap
           if options.chr:
               snps=[ x for x in getMap(options.chr,xleft=options.posleft,xright=options.posright)]
           else:
               snps=myMap.input_order[:]
       else:
           myMap=None
           snps=None
       dataset=parsePedFile(options.ped,params=params,snps=snps)
       if myMap is None:
           ## create a dummy map
           myMap=data.Map()
           for i,mysnp in enumerate(sorted(dataset.snp.keys(),key = lambda x:dataset.snpIdx[x])):
               myMap.addMarker(M=mysnp,C='0',posG=i,posP=i)
       return {'dataset':dataset,'map':myMap}
    ## 2. --file : matching ped and map files
    if options.file:
        mapFic=options.file+'.map'
        pedFic=options.file+'.ped'
        ## get Map and SNP info
        myMap=parseMapFile(mapFic)
        if options.chr:
            mySnpIdx=_get_snpIdx(myMap,options.chr,options.posleft,options.posright,options.other_map)
        # if options.other_map:
        #     getMap=myMap.genetMap
        # else:
        #     getMap=myMap.physMap
        # if options.chr:
        #     mySnps=[ x for x in getMap(options.chr,xleft=options.posleft,xright=options.posright)]
        #     mySnpIdx=[]
        #     for i,s in enumerate(myMap.input_order):
        #         if s in mySnps:
        #             mySnpIdx.append(i)
        else:
            mySnpIdx=None
        ## read Ped Data
        dataset=parsePedFile_nohead(pedFic,snpNames=myMap.input_order,params=params,snpidx=mySnpIdx)
        return {'dataset':dataset,'map':myMap}
    ## 3. --bfile : binary plink files
    if options.bfile:
        return parsePlinkBfile(options.bfile,params=params,options=options)
               
######################## PED / MAP FILES ##############################

def parsePedFile_nohead(fileName,snpNames,params=defaultParams,snpidx=None,excludeInd=[],noPheno=False):
    '''
    Read plink PED file
    
    Parameters
    --------------
    filename : PED file name
    snpNames : name of the SNPs in file (order matching)
    params : parameters for coding used in input files
    isnps : indices of SNPs to read, default is to read all SNPs
    excludeInd : list of IDs of individuals to exclude.
    noPheno : indicates if there is no phenotype column

    Returns
    ----------
    dataset : an instance of the data.Dataset class, 
              missing SNP information (needs to be added after call)

    See also
    ----------
    parsePedFile : for version with header line 
    data.Dataset : for help on the returned object
    '''
    print 'Reading',fileName
    ## count number of individuals
    with open(fileName) as f:
        for nindiv,l in enumerate(f):
            pass
    nindiv+=1
    ## manage SNPs
    nsnp=len(snpNames)
    if snpidx is None:
        snpidx=range(nsnp)
    ## create empty dataset
    fic=open(fileName)    
    dataset=data.Dataset(fileName,nsnp=len(snpidx),nindiv=nindiv)
    ## add SNP info
    for name in [snpNames[i] for i in snpidx]:
        dataset.addSnp(name)
    ## add individual info
    for ligne in fic:
        buf=ligne.split()
        sys.stdout.write('\tIND %16s\r'%buf[1])
        sys.stdout.flush()
        if buf[1] in excludeInd:
            continue
        if noPheno:
            phe=None
            genotype=buf[5:]
        else:
            phe=((buf[5]!=params['MissPheno']) and float(buf[5])) or None
            genotype=buf[6:]
    
        all1=genotype[::2]
        all2=genotype[1::2]
        try:
            assert len(all1)==len(snpNames) and len(all2)==len(snpNames)
        except AssertionError:
            print 'Genotype problem with ',buf[:6],len(all1),len(all2),len(snpNames)
            raise
        dataset.addIndividual(pop=buf[0],
                              ID=buf[1],
                              fatherID=((buf[2]!=params['MissParent']) and buf[2]) or None,
                              motherID=((buf[3]!=params['MissParent']) and buf[3]) or None,
                              sex=(buf[4]!=params['MissSex'] and buf[4]) or None,
                              phenotype=phe
                              )
        for s in snpidx:
            if all1[s]==params['MissGeno'] or all2[s]==params['MissGeno']:
                dataset.SetGenotype(buf[1],snpNames[s],missing)
            else:
                dataset.SetGenotype(buf[1],snpNames[s],all1[s]+all2[s])
        sys.stdout.flush()
    print
    return dataset

def parsePedFile_head(fileName,params=defaultParams,snps=None,excludeInd=[],noPheno=False):
    '''
    Read a genotype file (ped format with header line # [list of SNP names])
    options are:
      - snps: ordered list of snps for which to store information
      - excludeInd: don't store info for these individuals
    '''
    print 'Reading',fileName
    ## count number of individuals ( = # lines -1)
    with open(fileName) as f:
        for nindiv,l in enumerate(f):
            pass
    fic=open(fileName)    
    firstline=fic.readline().split()
    if firstline[0]!='#':
        print '\n\nFirst line has wrong format'
        print 'expected : # SNP1 SNP2 SNP3 ...'
        print 'got : ',' '.join(firstline[:5]),'...\n\n'
        raise ValueError
    ## construct a Hash of SNP names and positions
    allSNP={}
    isnp=0
    for snp in firstline[1:]:
        allSNP[snp]=isnp
        allSNP[key+str(isnp)]=snp
        isnp+=1
    ## now get indices of the requested snps in a sorted list
    nNotFound=0
    if snps is not None:
        snpidx=[]
        for SNPname in snps:
            try:
                insort(snpidx,allSNP[SNPname])
            except KeyError:
                ##print SNPname,'not found in file'
                nNotFound+=1
                continue
        print '# SNPs not found',nNotFound
    else:
        snpidx=range(len(firstline[1:]))
    dataset=data.Dataset(fileName,nsnp=len(snpidx),nindiv=nindiv)
    if snps is not None:
        ## Add SNPs in the order of snps argument
        for SNPname in snps:
            if not allSNP.has_key(SNPname):
                continue
            dataset.addSnp(SNPname)
    else:
        ## Add SNPs in the order of the file
        for SNPname in firstline[1:]:
            dataset.addSnp(SNPname)

    for ligne in fic:
        buf=ligne.split()
        sys.stdout.write('\tIND %16s\r'%buf[1])
        sys.stdout.flush()
        if buf[1] in excludeInd:
            continue
        if noPheno:
            phe=None
            genotype=buf[5:]
        else:
            phe=((buf[5]!=params['MissPheno']) and float(buf[5])) or None
            genotype=buf[6:]
    
        all1=genotype[::2]
        all2=genotype[1::2]
        try:
            assert len(all1)==len(allSNP)/2 and len(all2)==len(allSNP)/2
        except AssertionError:
            print 'Genotype problem with ',buf[:6],len(all1),len(all2),len(allSNP)/2
            continue
        dataset.addIndividual(pop=buf[0],
                              ID=buf[1],
                              fatherID=((buf[2]!=params['MissParent']) and buf[2]) or None,
                              motherID=((buf[3]!=params['MissParent']) and buf[3]) or None,
                              sex=(buf[4]!=params['MissSex'] and buf[4]) or None,
                              phenotype=phe
                              )
        for s in snpidx:
            if all1[s]==params['MissGeno'] or all2[s]==params['MissGeno']:
                dataset.SetGenotype(buf[1],allSNP[key+str(s)],missing)
            else:
                dataset.SetGenotype(buf[1],allSNP[key+str(s)],all1[s]+all2[s])
        sys.stdout.flush()
    print
    return dataset

parsePedFile=parsePedFile_head

def parsePedFile_fast(fileName,params=defaultParams,snps=None,excludeInd=[],noPheno=False):
    '''
    Read a genotype file (ped format with header line # [list of SNP names])
    options are:
      - snps: ordered list of snps for which to store information
      - excludeInd: don't store info for these individuals
    '''
    print 'Reading',fileName
    ## count number of individuals (=  number of lines -1)
    with open(fileName) as f:
        for nindiv,l in enumerate(f):
            pass
    fic=open(fileName)    
    firstline=fic.readline().split()
    if firstline[0]!='#':
        print '\n\nFirst line has wrong format'
        print 'expected : # SNP1 SNP2 SNP3 ...'
        print 'got : ',' '.join(firstline[:5]),'...\n\n'
        raise ValueError
    ## construct a Hash of SNP names and positions
    allSNP={}
    isnp=0
    for snp in firstline[1:]:
        allSNP[snp]=isnp
        allSNP[key+str(isnp)]=snp
        isnp+=1
    ## now get indices of the requested snps in a sorted list
    nNotFound=0
    if snps is not None:
        snpidx=[]
        for SNPname in snps:
            try:
                insort(snpidx,allSNP[SNPname])
            except KeyError:
                ##print SNPname,'not found in file'
                nNotFound+=1
                continue
        print '# SNPs not found',nNotFound
    else:
        snpidx=range(len(firstline[1:]))
    dataset=data.Dataset(fileName,nsnp=len(snpidx),nindiv=nindiv)
    if snps is not None:
        ## Add SNPs in the order of snps argument
        for SNPname in snps:
            if not allSNP.has_key(SNPname):
                continue
            dataset.addSnp(SNPname)
    else:
        ## Add SNPs in the order of the file
        for SNPname in firstline[1:]:
            dataset.addSnp(SNPname)
    for ligne in fic:
        if noPheno:
            phe=None
            buf=ligne.split(None,5)
        else:
            buf=ligne.split(None,6)
            phe=((buf[5]!=params['MissPheno']) and float(buf[5])) or None
        #print buf[0:10]
        sys.stdout.write('\tIND %16s\r'%buf[1])
        sys.stdout.flush()
        if buf[1] in excludeInd:
            continue
        genotype=re.sub(r'\s','',buf[-1])
        all1=genotype[::2]
        all2=genotype[1::2]
        try:
            assert len(all1)==len(allSNP)/2 and len(all2)==len(allSNP)/2
        except AssertionError:
            print 'Genotype problem with ',buf[:6],len(all1),len(all2),len(allSNP)/2
            continue
        dataset.addIndividual(pop=buf[0],
                              ID=buf[1],
                              fatherID=((buf[2]!=params['MissParent']) and buf[2]) or None,
                              motherID=((buf[3]!=params['MissParent']) and buf[3]) or None,
                              sex=(buf[4]!=params['MissSex'] and buf[4]) or None,
                              phenotype=phe
                              )
        for s in snpidx:
            if all1[s]==params['MissGeno'] or all2[s]==params['MissGeno']:
                dataset.SetGenotype(buf[1],allSNP[key+str(s)],missing)
            else:
                dataset.SetGenotype(buf[1],allSNP[key+str(s)],all1[s]+all2[s])
        sys.stdout.flush()
    print
    return dataset
   
def parseMapFile(fileName):
    fic=open(fileName)
    print 'Reading',fileName
    myMap=data.Map()
    for ligne in fic:
        buf=ligne.split()
        try:
            mychr=int(buf[0])
        except:
            mychr=buf[0]
        mychr=buf[0]
        myMap.addMarker(M=buf[1],C=mychr,posG=float(buf[2]),posP=float(buf[3]))
    return myMap

################ BED / BIM / FAM files #############################

def parsePlinkBfile(prefix,noPheno=False,params=defaultParams,options=None):
    fam=prefix+'.fam'
    bim=prefix+'.bim'
    bed=prefix+'.bed'
    ## individual data
    idata=parseFamFile(fam) ## options not implemented
    ni=len(idata)
    ### SNP data
    sdata=parseBimFile(bim)
    if options and options.chr:
        mySnpIdx=_get_snpIdx(sdata['map'],options.chr,options.posleft,options.posright,options.other_map)
        ns=len(mySnpIdx)
    else:
        ns=len(sdata['snps'])
        mySnpIdx=range(ns)
    dataset=data.Dataset(prefix,nsnp=ns,nindiv=ni)
    ## fill in SNP data
    for s in [sdata['snps'][i] for i in mySnpIdx]:
        dataset.addSnp(s.name)
        dataset.snp[s.name].initAlleles(s.alleles[0],s.alleles[1])
    ## fill in indiv data
    for ind in idata:
        dataset.addIndividual(pop=ind[0],
                              ID=ind[1],
                              fatherID=ind[2],
                              motherID=ind[3],
                              sex=ind[4],
                              phenotype=ind[5])
    fillBedData(bed,dataset.Data,mySnpIdx)
    return {'dataset':dataset,'map':sdata['map']}

def parseFamFile(fileName,noPheno=False,params=defaultParams):
    indiv_data=[]
    with open(fileName) as f:
        for ligne in f:
            buf=ligne.split()
            if noPheno:
                phe=None
            else:
                phe=((buf[5]!=params['MissPheno']) and float(buf[5])) or None
  
            indiv_data.append([buf[0], \
                               buf[1],  \
                               ((buf[2]!=params['MissParent']) and buf[2]) or None, \
                               ((buf[3]!=params['MissParent']) and buf[3]) or None, \
                               (buf[4]!=params['MissSex'] and buf[4]) or None, \
                               buf[5]])
    return indiv_data

def parseBimFile(fileName):
    SNPs=[]
    myMap=data.Map()
    with open(fileName) as f:
        for ligne in f:
            buf=ligne.split()
            ## create SNP
            myS=data.SNP(buf[1])
            myS.alleles[0]=buf[4]
            myS.alleles[1]=buf[5]
            ## add marker to map
            try:
                mychr=int(buf[0])
            except:
                mychr=buf[0]
            mychr=buf[0]
            myMap.addMarker(M=buf[1],C=mychr,posG=float(buf[2]),posP=float(buf[3]))
            SNPs.append(myS)
    return {'snps':SNPs,'map':myMap}


def fillBedData(fileName,DataMatrix,snpidx=None):
    n_indiv,n_snp=DataMatrix.shape
    ## use dict type for faster look up
    ## key index in total, val : index in DataMatrix
    if snpidx is None:
        snpidx=dict([ (i,i) for i in range(n_snp)])
    else:
        snpidx=dict([ (s,i) for i,s in enumerate(snpidx)])
    cur_snp=0
    ## bit enumerator
    def bits(f):
        bytes = (ord(b) for b in f.read())
        for b in bytes:
            for i in xrange(8):
                yield (b >> i) & 1
    ## bit pair to geno converter
    def bpair_2_geno(pair):
        '''
        from plink doc
        '''
        if pair=='00':
            return 0
        if pair=='11':
            return 2
        if pair=='01':
            return 1
        if pair=='10':
            return missing
    ## constants
    target_magic='0011011011011000'
    snp_major='10000000'    
    indiv_major='00000000'
    ### let's go
    bed_stream=open(fileName, 'r')
    magic=[]
    f_mode=[]
    i_major=0
    i_minor=0
    proceed_2_next_bit=False
    bit_geno=[]
    nbit=-1
    for b in bits(bed_stream):
        nbit+=1
        ## Test for magic number
        if nbit<16:
            magic.append(str(b))
            continue
        if nbit==16:
            magic=''.join(magic)
            if magic!=target_magic:
                print 'Not a bed file'
                raise ValueError
        ## Find out packing mode
        if nbit<24:
            f_mode.append(str(b))
            continue
        if nbit==24:
            f_mode=''.join(f_mode)
            if f_mode==snp_major:
                pack_len=n_indiv 
            elif f_mode==indiv_major:
                pack_len=n_snp
            else:
                print 'Cannot determine packing mode, abort'
                raise ValueError
        ## Now i>24, we are reading genotypes
        ## case where we just wait for next bit to start reading again
        if proceed_2_next_bit:
            if nbit%8!=0:
                continue
            else:
                proceed_2_next_bit=False
        ## reading genotypes        
        bit_geno.append(str(b))
        if nbit%2!=0:
            ##assert len(bit_geno)==2
            ## print nbit,i_major,i_minor,bit_geno
            geno=bpair_2_geno(''.join(bit_geno))
            if f_mode==snp_major:
                try:
                    data_idx=snpidx[i_major]
                    DataMatrix[i_minor,data_idx]=geno
                except KeyError:
                    pass
            else:
                try:
                    data_idx=snpidx[i_minor]
                    DataMatrix[i_major,data_idx]=geno
                except KeyError:
                    pass
            i_minor += 1 
            bit_geno=[]
        if i_minor==pack_len:
            proceed_2_next_bit=True
            i_minor=0
            i_major+=1
        
