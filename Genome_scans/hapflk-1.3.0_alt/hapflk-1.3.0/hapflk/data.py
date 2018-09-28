from __future__ import division
from operator import itemgetter
from bisect import *
import sys
import numpy as np
##import numpy.ma as ma
from hapflk import missing, complete_cases, is_na

zero=np.int16(0)
un=np.int16(1)
deux=np.int16(2)

class Dataset():
    def __init__(self,fileName,nsnp,nindiv):
        self.name=fileName
        ## Individuals
        self.nindiv=nindiv
        self.iindiv=0
        self.indiv={}
        ## Individual index in Data Matrix
        self.indivIdx={}
        ## SNPs
        self.nsnp=nsnp
        self.isnp=0
        self.snp={}
        ## SNP indices in Data Matrix
        self.snpIdx={}
        self._initDataMatrix()
        ## populations
        ## populations are accessed by a name (key of the dict)
        ## the value is a vector of booleans taking a True value for individuals
        ## belonging to the pop otherwise
        ## These vectors can be used for fast access to population genotypes e.g.:
        ## vpop=self.populations['mypop']
        ## mypop_genotypes=self.Data[vpop,]
        self.populations={}
    def addSnp(self,Name):
        ''' Add a snp to the dataset '''
        if self.snp.has_key(Name):
            return
        self.snp[Name]=SNP(Name)
        self.snpIdx[Name]=self.isnp
        self.isnp+=1
    def rmSnp(self,Name):
        ''' 
        Remove a snp from the dataset 
        The SNP is made unaccessible but the data stay in the Data Matrix
        '''
        if not self.snp.has_key(Name):
            return
        else:
            del self.snp[Name]
            del self.snpIdx[Name]
            del self._allSnps
    def addIndividual(self,ID,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        ''' Add an individual to the dataset'''
        if self.indiv.has_key(ID):
            print "Warning: Individual ",ID," has been seen before !"
            return
        self.indiv[ID]=Individual(ID,pop,sex,fatherID,motherID,phenotype)
        self.indivIdx[ID]=self.iindiv
        if pop is not None:
            try:
                vpop=self.populations[pop]
            except KeyError:
                vpop=np.array([0]*self.nindiv,dtype=bool)
                self.populations[pop]=vpop
            vpop[self.iindiv]=True
        self.iindiv+=1
    def _initDataMatrix(self):
        ''' Initialize Data Matrix to an int16 array filled with -1'''
        self.Data=-np.ones(shape=(self.nindiv,self.nsnp),dtype='int16')
    def IndivGenotype(self,indiv,snps=None):
        '''
        Return (0,1,2) genotype of individual INDIV at all SNPs in the dataset
        '''
        try:
            i=self.indivIdx[indiv]
        except KeyError:
            print 'Individual',indiv,'not in Data'
            raise
        if snps is None:
            return [ (x<0 and missing) or x for x in  self.Data[i,]]
        else:
            return [ (x<0 and missing) or x for x in self.Data[i,[self.snpIdx[s] for s in snps]]]
    def IndivGenotype_char(self,indiv,snps=None):
        igeno=self.IndivGenotype(indiv,snps)
        cgeno=[]
        if snps is None:
            snps=self.AllSnps()
        for i,s in enumerate(snps):
            cgeno.append(self.Alleles(igeno[i],s))
        return cgeno
    def AllSnps(self):
        try:
            return self._allSnps
        except AttributeError:
            self._allSnps=sorted(self.snp,key=lambda x: self.snpIdx[x])
        return self._allSnps
    def SnpGenotype(self,snp):
        '''
        Return genotype of snp SNP for all individuals in the dataset
        '''
        try:
            j=self.snpIdx[snp]
        except KeyError:
            return [missing]*self.nindiv
        return [ (x<0 and missing) or x for x in self.Data[:,j]]
    def Alleles(self,geno,snp):
        '''
        Return the alleles corresponding to genotype GENO at snp SNP
        '''
        if geno is missing:
            return ('?','?')
        snpAll=self.snp[snp].alleles
        if geno==0:
            return (snpAll[0],snpAll[0])
        elif geno==1:
            return (snpAll[0],snpAll[1])
        else:
            assert geno==2
            return (snpAll[1],snpAll[1])
    def Genotype(self,indiv,snp):
        '''
        Return the genotype of individual INDIV  at snp SNP 
        If SNP is not found in the dataset, returns numpy.ma.masked
        '''
        try:
            i=self.indivIdx[indiv]
        except KeyError:
            print 'Individual',indiv,'not in Data'
            raise
        try:
            j=self.snpIdx[snp]
        except KeyError:
            return missing
        return self.Data[i,j]<0 and missing or self.Data[i,j]
    def SetGenotypeSlow(self,indiv,snp,value):
        '''
        Set the genotype of individual INDIV at snp SNP to value VALUE
        Correct Genotype is a value of length 2, with alleles matching the SNP
        '''
        try:
            i=self.indivIdx[indiv]
        except KeyError:
            print 'Individual',indiv,'not in Data'
            raise
        try:
            j=self.snpIdx[snp]
        except KeyError:
            print 'SNP',snp,'not in Data'
            raise
        if self.Data[i,j] > -1:
            self.snp[snp].delGenotype(self.Data[i,j])
        geno=self.snp[snp].addObservation(value)
        if geno is not missing:
            self.indiv[indiv].callrate+=1
        self.Data[i,j]=geno
    def SetGenotype(self,indiv,snp,value):
        '''
        Fast Version of SetGenotype with no checking for errors 
        Use with caution
        --
        Set the genotype of individual INDIV at snp SNP to value VALUE
        Correct Genotype is a value of length 2, with alleles matching the SNP
        '''
        try:
            i=self.indivIdx[indiv]
            j=self.snpIdx[snp]
            self.Data[i,j]=self.snp[snp].addObservation(value)     
        except:
            print 'Error'
            raise
    def printSNP(self,stream=sys.stdout):
        '''
        Print SNP information
        '''
        print >>stream,'name\tA1\tA2\tcallRate'
        for sName in sorted(self.snp,key=lambda x:self.snpIdx[x]):
            geno=self.SnpGenotype(sName)
            s=self.snp[sName]
            print >>stream,s.name,s.alleles[0],s.alleles[1],np.mean([g is not missing for g in geno])
    ##### PEDIGREE AND FAMILIES
    def BuildPedigree(self):
        '''
        Construct the pedigree relationships between individuals in the dataset
        Identify founders and parent-offspring links
        '''
        self.pedigree=Pedigree(self.indiv.values())
        self.founders=np.array([0]*self.nindiv,dtype=bool)
        for guy in self.indiv.values():
            founder=True
            if self.pedigree.setFather(guy,guy.fatherID):
                founder=False
            if self.pedigree.setMother(guy,guy.motherID):
                founder=False
            if founder:
                self.pedigree.founders.append(guy)
                self.founders[self.indivIdx[guy.id]]=True
    def BuildFamilies(self):
        '''
        Cluster related individuals in families
        Sets the self.families attribute to a list of the families
        '''
        try:
            ped=self.pedigree
        except AttributeError:
            self.BuildPedigree()
            ped=self.pedigree
        ped.BuildFamilies()
        self.families=ped.families
    def UnrelatedIndividuals(self):
        try:
            ped=self.pedigree
        except AttributeError:
            self.BuildPedigree()
            ped=self.pedigree
        return ped.UnrelatedIndividuals()
    def NuclearFamilies(self):
        '''
        Generator that yields all nuclear families in the dataset (with both parents present)
        '''
        try:
            allFams=self.families
        except AttributeError:
            self.BuildFamilies()
            allFams=self.families
        for fam in allFams:
            for nucFam in fam.NuclearFamilies():
                yield nucFam
    
    def HalfSibFamilies(self):
        '''
        Generator that yields all half-sib families in the dataset.
        '''
        try:
            allFams=self.families
        except AttributeError:
            self.BuildFamilies()
            allFams=self.families
        for fam in allFams:
            for hsFam in fam.HalfSibFamilies():
                yield hsFam
    
    def checkSegregation(self,autosome=True,sex_type='XY'):
        # if not autosome:
        #     if sex_type=='XY':
        #         return self.checkSegregationXY()
        #     elif sex_type=='ZW':
        #         return self.checkSegregationZW()
        try:
            allFams=self.families
        except AttributeError:
            self.BuildFamilies()
            allFams=self.families
        self.MendelErrors={}
        for fam in allFams:
            for node in fam.nonFounders:
                o=node.indiv
                ## father <-> offspring
                if (autosome or sex_type=='ZW') and node.father is not None:
                    father=node.father.indiv
                    for snp in self.AllSnps():
                        go=self.Genotype(o.id,snp)
                        gf=self.Genotype(father.id,snp)
                        if go is missing or gf is missing:
                            continue
                        if (go==0 and gf==2) or (go==2 and gf==0):
                            try:
                                self.SetGenotype(o.id,snp,missing)
                                self.SetGenotype(father.id,snp,missing)
                            except:
                                print o.id,father.id,go,gf
                                
                            try :
                                self.MendelErrors[snp]+=1
                            except KeyError:
                                self.MendelErrors[snp]=1
                            try :
                                self.MendelErrors[o.id]+=1
                            except KeyError:
                                self.MendelErrors[o.id]=1
                            try :
                                self.MendelErrors[father.id]+=1
                            except KeyError:
                                self.MendelErrors[father.id]=1
                ## mother <-> offspring
                if (autosome or sex_type=='XY') and node.mother is not None:
                    mother=node.mother.indiv
                    for snp in self.AllSnps():
                        go=self.Genotype(o.id,snp)
                        gm=self.Genotype(mother.id,snp)
                        if go is missing or gm is missing:
                            continue
                        if (go==0 and gm==2) or (go==2 and gm==0):
                            self.SetGenotype(o.id,snp,missing)
                            self.SetGenotype(mother.id,snp,missing)
                            try :
                                self.MendelErrors[snp]+=1
                            except KeyError:
                                self.MendelErrors[snp]=1
                            try :
                                self.MendelErrors[o.id]+=1
                            except KeyError:
                                self.MendelErrors[o.id]=1
                            try :
                                self.MendelErrors[mother.id]+=1
                            except KeyError:
                                self.MendelErrors[mother.id]=1
        return self.MendelErrors
    ### POPULATION GENETICS methods
    def compute_pop_frq(self,verbose=True):
        '''
        Compute allele frequencies for each population in the dataset, considering only
         founders.
        '''
        try:
            fvec=self.founders
        except AttributeError:
            self.BuildPedigree()
            fvec=self.founders
        pop_names=[]
        frqs=[]
        for pname,pvec in self.populations.items():
            if verbose:
                sys.stdout.write('\t %16s\r'%pname)
                sys.stdout.flush()
            pop_founders=pvec&fvec
            w=complete_cases(self.Data[pop_founders,])
            try:
                frqs.append(0.5*np.average(self.Data[pop_founders,],axis=0,weights=w))
            except ZeroDivisionError:
                frqs.append(0.5*np.ma.average(self.Data[pop_founders,],axis=0,weights=w))
            pop_names.append(pname)
        print
        return {"pops":pop_names,'freqs':np.vstack(frqs)}
   
class Family():
    def __init__(self,founders=[],nonFounders=[]):
        self.founders=founders
        self.nonFounders=nonFounders
    def __str__(self):
        tw='Founders: '+'x'.join([x.indiv.id for x in self.founders])+' with '+str(len(self.nonFounders))+' offspring'
        return tw
    def addFounder(self,founder):
        ''' Add a node to the founders '''
        if founder not in self.founders:
            self.founders.append(founder)
    def addNonFounders(self,offspring):
        ''' Add a list of nodes to the offspring '''
        if offspring not in self.nonFounders:
            self.nonFounders.append(offspring)
    def NuclearFamilies(self):
        '''
            Build Nuclear Families within an extended pedigree (family)
            Family founders are ordered so that fam.founders[0] is dad
        '''
        FSfam=[]
        children=list(self.nonFounders)
        while len(children)>0:
            firstFS=children[0]
            papa=firstFS.father
            maman=firstFS.mother
            if papa is not None and maman is not None:
                enfants=[off for off in children if off.father==papa and off.mother==maman]
                FSfam.append(Family(founders=[papa,maman],nonFounders=enfants))
                for off in enfants:
                    children.remove(off)
            else:
                children.remove(firstFS)
        return FSfam
    def _rmHiddenFullSibs(self,enfants,fatherConnected=True):
        '''
        remove "hidden" full sibs from a set of individuals (i.e.
        have same motherID (fatherConnected=True) 
        or fatherID (fatherConnected=False))
        '''
        keepThem=[]
        trashThem=[]
        otherParents={}
        ## make sure we consider enfants always in the same order
        enfants=sorted(enfants,key=lambda x:x.indiv.id)
        for node in sorted(enfants,key=lambda x:x.indiv.callrate,reverse=True):
            guy=node.indiv
            if fatherConnected:
                otherID=guy.motherID
            else:
                otherID=guy.fatherID
            try:
                # we have encountered a full sib before
                # we don't keep that guy
                otherParents[otherID]+=1
                trashThem.append(node)
            except KeyError:
                otherParents[otherID]=1
                keepThem.append(node)
        return keepThem,trashThem
                    
    def HalfSibFamilies(self):
        '''
        Return HalfSib Families contained in the family
        '''
        HSfam=[]
        children=list(self.nonFounders)
        while len(children)>0:
            firstHS=children[0]
            papa=firstHS.father
            maman=firstHS.mother
            if papa is not None:
                if maman is None:
                    ## be sure all of them have none mothers
                    enfants=[off for off in children if off.father==papa and off.mother is None]
                    for off in enfants:
                        children.remove(off)
                    ## remove hidden FS
                    ##enfants,rmvd=self._rmHiddenFullSibs(enfants)
                    HSfam.append(Family(founders=[papa],nonFounders=enfants))
                else:
                    ## nuclear family: remove everyone
                    enfants=[off for off in children if off.father==papa and off.mother==maman]
                    for off in enfants:
                        children.remove(off)
            elif maman is not None:
                ## father is None
                enfants=[off for off in children if off.mother==maman and off.father is None]
                for off in enfants:
                    children.remove(off)
                ## remove hidden FS
                ##enfants,rmvd=self._rmHiddenFullSibs(enfants,fatherConnected=False)
                HSfam.append(Family(founders=[maman],nonFounders=enfants))
            else:
                children.remove(firstHS)
        return HSfam
    
    
class Pedigree():
    ''' Pedigree relationships between individuals'''
    def __init__(self,guys=[]):
        self.nodes={}
        self.founders=[]
        self.families=[]
        for guy in guys:
            self.addNode(guy)
    def addNode(self,guy):
        ''' Add a new guy to the pedigree'''
        self.nodes[guy.id]=PedNode(guy)
    def setFather(self,son,dadID):
        ''' 
        Set that dadId is dad of sonId 
        return True if dad found
        return False if not found
        raise if sonID not found
        '''
        if not self.nodes.has_key(son.id):
            print "Don't know that guy",son.id
            raise KeyError 
        try:
            self.nodes[son.id].father=self.nodes[dadID]
            self.nodes[dadID].children.append(self.nodes[son.id])
            return True
        except KeyError:
            return False
    def setMother(self,son,momID):
        ''' 
        Set that momId is mom of sonId 
        return True if mom found
        return False if not found
        raise if sonID not found
        '''
        if not self.nodes.has_key(son.id):
            print "Don't know that guy",sonId
            raise KeyError 
        try:
            self.nodes[son.id].mother=self.nodes[momID]
            self.nodes[momID].children.append(self.nodes[son.id])
            return True
        except KeyError:
            return False    
    def _getRelatives(self,node,connected):
        if node in connected:
            return
        connected.append(node)
        for off in [x for x in node.children]:
            self._getRelatives(off,connected)
        if node.father!=None:
            self._getRelatives(node.father,connected)
        if node.mother!=None and node.mother not in connected:
            self._getRelatives(node.mother,connected)
    def UnrelatedIndividuals(self):
        return [x.id for x in self.founders if len(self.nodes[x.id].children)==0]
    def BuildFamilies(self):
        '''
        Cluster individuals by families
        '''
        founders2go=[self.nodes[x.id] for x in self.founders[:]]
        offspring2go=list(set(self.nodes.values())-(set(founders2go)&set(self.nodes.values())))
        while len(founders2go)>0:
            relatives=[]
            self._getRelatives(founders2go[0],relatives)
            famFounders=[]
            famNonFounders=[]
            for x in relatives:
                if x in founders2go:
                    famFounders.append(x)
                else:
                    famNonFounders.append(x)
            if len(famNonFounders)>0:
                self.families.append(Family(founders=famFounders,nonFounders=famNonFounders))
            for x in famFounders:
                founders2go.remove(x)
            for x in famNonFounders:
                try:
                    offspring2go.remove(x)
                except ValueError:
                    print 'pula',x.indiv.id

class PedNode():
    def __init__(self,indiv,father=None,mother=None):
        self.indiv=indiv
        self.father=father
        self.mother=mother
        self.children=[]
    def addChild(self,child):
        if child not in self.children:
            self.children.append(child)
        
class Individual():
    def __init__(self,id,pop=None,sex=None,fatherID=None,motherID=None,phenotype=None):
        self.id=id
        self.pop=pop
        self.sex=sex
        self.fatherID=fatherID
        self.motherID=motherID
        self.phenotype=[phenotype]
        self.callrate=0


class SNP():
    def __init__(self,Name):
        self.name=Name
        self.alleles=[None,None] ## maps 0,1 to alleles
        self.mapping={} ## maps allele to 0,1
        self.counts=np.zeros(3,dtype='int16') ## number(a1),number(a2),number(missing)
        self.Gcounts=[0,0,0]
        self.addObservation=self.addObservationNotFull
    def callrate(self):
        return 0.5*(self.counts[0]+self.counts[1])/(0.5*(self.counts[0]+self.counts[1])+self.counts[2])
    def initAlleles(self,a1,a2):
        self.alleles[0]=a1
        self.alleles[1]=a2
        self.mapping[a1]=0
        self.mapping[a2]=1
        self._create_mapping()
        self.addObservation=self.addObservationFull
    def delGenotype(self,geno):
        '''
        Remove a genotype from observed alleles
        '''
        if geno == missing:
            self.counts[2] -= 1
        else:
            self.counts[0] -= 2-geno
            self.counts[1] -= geno
    def addGenotype(self,geno):
        '''
        Add a genotype to observed alleles
        '''
        if geno==missing:
            self.counts[2]-=1
        else:
            self.counts[0] += 2-geno
            self.counts[1] += geno
    def _create_mapping(self):
        self.g2int={}
        self.g2int[self.alleles[0]+self.alleles[0]]=zero
        self.g2int[self.alleles[0]+self.alleles[1]]=un
        self.g2int[self.alleles[1]+self.alleles[0]]=un
        self.g2int[self.alleles[1]+self.alleles[1]]=deux
    def addObservationFull(self,value):
        try:
            return self.g2int[value]
        except KeyError:
            return missing
    def addObservationNotFull(self,value):
        ''' 
        Add an observed genotype ('AA','01'...) to the SNP 
        Return the 'integer genotype',0,1,2
        if two alleles have already been observed, switch to fast
        version
        '''
        geno=missing
        try:
            i1,i2=value[0],value[1]
        except :
            return geno
        try:
            g1=self.mapping[i1]
        except KeyError:
            if self.alleles[0] is None:
                self.alleles[0]=i1
                self.mapping[i1]=zero
                g1=zero
            elif self.alleles[1] is None:
                self.alleles[1]=i1
                self.mapping[i1]=un
                g1=un
                self._create_mapping()
                self.addObservation=self.addObservationFull
            else:
                return geno
        try:
            g2=self.mapping[i2]
        except KeyError:
            assert self.alleles[0]!=None
            if self.alleles[1] is None:
                self.alleles[1]=i2
                self.mapping[i2]=un
                g2=un
                self._create_mapping()
                self.addObservation=self.addObservationFull
            else:
                return geno
        geno=g1+g2
        return geno

class Locus():
    def __init__(self,nom,position):
        self.name=nom
        self.pos=position
    def __cmp__(self,other):
        return cmp(self.pos,other.pos)

class Map():
    def __init__(self):
        '''
        Map class contains maps (duh)
        The chromosome names are stored in the chromosomes attribute (in sorted order)
        Each chromosome has two maps (genetic/RH and physical)
        '''
        self.chromosomes=[] 
        self.cartes={}
        self.markers={}
        self.input_order=[]
    def write(self,stream=sys.stdout):
        for chrom in self.chromosomes:
            for mk in self.physMap(chrom):
                pos=self.position(mk)
                print >>stream,'\t'.join([str(x) for x in [pos[0],mk,int(pos[1]),int(pos[2])]])
        return
    def addMarker(self,M,C,posG=0,posP=0):
        '''
        Add a marker on chromosome C, with genetic position posG and
        physical position posP
        '''
        self.input_order.append(M)
        if not self.cartes.has_key(C):
            insort(self.chromosomes,C)
            self.cartes[C]={'Genet':[],'Phys':[]}
        locG=Locus(M,posG)
        locM=Locus(M,posP)
        self.markers[M]=(C,locG,locM)
        insort(self.cartes[C]['Genet'],locG)
        insort(self.cartes[C]['Phys'],locM)
    def delMarker(self,M):
        try:
            locus = self.markers[M]
        except KeyError:
            print "No such marker"
        del self.markers[M]
        self.cartes[locus[0]]['Genet'].remove(locus[1])
        self.cartes[locus[0]]['Phys'].remove(locus[2])
    def genetMap(self,C,xleft=-1,xright=np.inf):
        '''
        Returns and iterator on the genetic map between position xleft and xright 
        on chromosome C
        '''
        try:
            carte=self.cartes[C]['Genet']
        except KeyError:
            return []
        mk=carte[bisect_left(carte,Locus(0,xleft)):bisect_right(carte,Locus(0,xright))]
        return (m.name for m in mk)
    def physMap(self,C,xleft=0,xright=np.inf):
        '''
        Returns and iterator on the genetic map between position xleft and xright 
        on chromosome C
        '''
        try:
            carte=self.cartes[C]['Phys']
        except KeyError:
            return []
        mk=carte[bisect_left(carte,Locus(0,xleft)):bisect_right(carte,Locus(0,xright))]
        return (m.name for m in mk)
    def length(self,C):
        '''
        Return the length of maps of chromosome C [physical,genetic]
        '''
        carte=self.cartes[C]['Phys']
        lenP=carte[-1].pos-carte[0].pos
        carte=self.cartes[C]['Genet']
        lenG=carte[-1].pos-carte[0].pos
        return [lenP,lenG]
    def position(self,M):
        '''
        Returns chromosome,genetic,physical position of marker M
        '''
        try:
            locus=self.markers[M]
            return locus[0],locus[1].pos,locus[2].pos
        except KeyError:
            print M,'not found'
            return None
    def sort_loci(self,loci,physical=True):
        '''
        Sort loci according to their position on the map

        Parameters:
        ---------------

        loci : a list of loci names
        physical : if True sort according to physical position, o.w. use genetic  / RH position

        Return Value:
        -----------------

        list of the same loci, sorted
        '''
        new_list=[]
        for locus in loci:
            try:
                tmp=self.markers[locus]
                if physical:
                    myloc=tmp[2]
                else:
                    myloc=tmp[1]
                try:
                    mychr=int(tmp[0])
                except:
                    mychr=tmp[0]
            except KeyError:
                continue
            new_list.append((locus,mychr,myloc.pos))
        new_list=sorted(new_list,key=itemgetter(1,2))
        return [x[0] for x in new_list]
    
    def iterate_pos_windows(self,wlen,padlen=0):
        for c in self.chromosomes:
            pos_left=0
            pos_right=-1
            clen=self.length(c)[0]
            print c,clen
            while pos_right < clen:
                pos_right=pos_left+wlen
                yield (c,(pos_left-padlen)>0 and pos_left-padlen or 0,(pos_right+padlen)>clen and clen+1 or pos_right+padlen )
                pos_left=pos_right
                
