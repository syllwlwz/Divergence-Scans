import sys
import os
import tempfile
import numpy as np
from scipy.stats import chi2,norm
import nj

debug=False
if not debug:
    np.seterr(all='ignore')

def reynolds_multi_allelic(Mf,nloc,decompose=False):
    '''
    Compute Reynolds distances from multiallelic loci

    Parameters:
    -----------
    
    Mf : numpy array of allele frequencies in populations
         rows are populations
         columns are allele frequencies, by loci
    nloc : number of loci

    Returns:
    --------

    dist : numpy array (n x n) of reynolds genetic distances
    numerator    : numerator of the Reynolds Distance
    denominator  : denominator of the Reynolds Distance
    '''
    npop=Mf.shape[0]
    A=np.dot(Mf,Mf.T)
    dist=np.zeros((npop,npop))
    if decompose:
        numerator=np.zeros((npop,npop))
        denominator=np.zeros((npop,npop))
    else:
        numerator=None
        denominator=None
    for i in range(npop-1):
        for j in range(i+1,npop):
            dist[i,j]=0.5*(A[i,i]+A[j,j]-2*A[i,j])/(nloc-A[i,j])
            if decompose:
                numerator[i,j]=A[i,i]+A[j,j]-2*A[i,j]
                denominator[i,j]=2*(nloc-A[i,j])
    dist=dist+dist.T
    numerator=numerator+numerator.T
    denominator=denominator+denominator.T
    return dist,numerator,denominator
   
def reynolds(Mf):
    '''
    Compute Reynolds distances from SNP allele frequencies

    Parameters:
    ----------------

    Mf : numpy array of allele frequencies in populations
          rows are populations (n), columns are markers (p).

    Returns:
    -----------

    dist : numpy array (n x n) of reynolds genetic distances
    '''
    npop,nloc=Mf.shape
    dist=np.zeros((npop,npop))
    A=np.dot(Mf,Mf.T)+np.dot((1-Mf),(1-Mf).T)
    for i in range(npop-1):
        for j in range(i+1,npop):
            dist[i,j]=0.5*(A[i,i]+A[j,j]-2*A[i,j])/(nloc-A[i,j])
    dist=dist+dist.T
    return dist

def heterozygosity(Mf):
    '''
    Compute heterozygosity from SNP allele frequencies
    Parameters:
    ----------------

    Mf : numpy array of allele frequencies in populations
          rows are populations (n), columns are markers (p).

    Returns:
    -----------

    hzy : numpy array (n x 1) of mean heterozygosities
    '''
    npop,nloc=Mf.shape
    hzy_func = lambda x: 4*x^3*(1-x)+6*x^2*(1-x)^2+4*x*(1-x)^3
    return np.average(hzy_func(Mf),axis=1)
    
Fij_Rscript='''
midpoint= function(tree){
    require(phangorn)
    dm = cophenetic(tree)
    tree = unroot(tree)
    rn = max(tree$edge)+1
    maxdm = max(dm)
    ind =  which(dm==maxdm,arr=TRUE)[1,]
    tmproot = Ancestors(tree, ind[1], "parent")
    tree = phangorn:::reroot(tree, tmproot)
    edge = tree$edge
    el = tree$edge.length
    children = tree$edge[,2]
    left = match(ind[1], children)
    tmp = Ancestors(tree, ind[2], "all")
    tmp= c(ind[2], tmp[-length(tmp)])
    right = match(tmp, children)
    if(el[left]>= (maxdm/2)){
         edge = rbind(edge, c(rn, ind[1]))
         edge[left,2] = rn
         el[left] = el[left] - (maxdm/2)
         el = c(el, maxdm/2)
    }
    else{
        sel = cumsum(el[right])
        i = which(sel>(maxdm/2))[1]
        edge = rbind(edge, c(rn, tmp[i]))
        edge[right[i],2] = rn
        eltmp =  sel[i] - (maxdm/2)
        el = c(el, el[right[i]] - eltmp)
        el[right[i]] = eltmp
    }
    tree$edge.length = el
    tree$edge=edge
    tree$Nnode  = tree$Nnode+1
    phangorn:::reorderPruning(phangorn:::reroot(tree, rn))
}

kinship=function(D,outgroup="",keep_outgroup=FALSE)
  {
    require(ape)
    D=2*as.matrix(D)
    npop=sqrt(length(D))
    PopTree=nj(D)
    if (outgroup=="") {
      PopTree=midpoint(PopTree)
      npop=npop+1
    } else {
      PopTree=root(PopTree,outgroup)
      if (!keep_outgroup) {
       PopTree=drop.tip(PopTree,outgroup)
      } else {
       npop=npop+1
      }
    }
    ## negative branch length are meaningless in our context
    PopTree$edge.length=abs(PopTree$edge.length)
    ## get the tree structure as : [ father, son, length]
    edges=cbind(PopTree$edge,PopTree$edge.length)
    #Identify ancestral node as father node with no father itself
    father.nodes=unique(edges[,1])
    son.nodes=unique(edges[,2])
    ancestral=father.nodes[which(is.na(match(father.nodes,son.nodes)))]
    ## Now compute F matrix
    Fij=matrix(0,nrow=npop-1,ncol=npop-1,
    dimnames=list(PopTree$tip.label,PopTree$tip.label))
    branch.length=dist.nodes(PopTree)
    tips=1:(npop-1)
    route=vector("list",npop-1)
    for (ipop in tips) {
      ## First Fii = distance to the root
      Fij[ipop,ipop]=branch.length[ipop,ancestral]
      ## build route to root
      father=ipop
      while (father!=ancestral) {
        route[[ipop]]=c(route[[ipop]],father)
        edj=which(edges[,2]==father)
        father=edges[edj,1]
      }
      route[[ipop]]=c(route[[ipop]],ancestral)
    }
    ## Now compute the Fij i!= j
    for (ipop in tips[-length(tips)]) {
      for (jpop in (ipop+1):(npop-1)) {
        common.route=intersect(route[[ipop]],route[[jpop]])
        if (length(common.route)>1) {
          for (i in 2:length(common.route)) {
            n1=common.route[i-1]
            n2=common.route[i]
            Fij[ipop,jpop] = Fij[ipop,jpop] + branch.length[n1,n2]
          }
        }
        Fij[jpop,ipop]=Fij[ipop,jpop]
      }
    }
    return(list(matrix=Fij,pops=colnames(Fij),tree=PopTree))
  }
'''

def popKinship_fromFile(file_name,popnames):
    '''
    Read population kinship matrix from file

    File format is :

    Pop1 a11 a12 ...
    Pop2 a21 a22 ....

    where PopN is the population name (must be in popnames)
    aij is the kinship coefficient between pop i and j

    Parameters:
    -----------
    file_name : name of the file to read kinship from
    popnames : name of the pop to find in the file (order is preserved)

    See Also:
    ---------
    popKinship for estimating the kinship matrix
    '''

    data=[x.split() for x in open(file_name).readlines()]
    ##popnames2=popnames[:]
    # if outgroup is not None:
    #     popnames2.remove(outgroup)
    ordre=[popnames.index(x[0]) for x in data]
    if len(ordre) != len(data):
        print 'Not Enough populations in kinship file !'
        raise ValueError
    FF=np.zeros((len(popnames),len(popnames)),dtype=float)
    for i in range(len(popnames)):
        ix=ordre[i]
        for j in range(len(popnames)):
            jx=ordre[j]
            FF[ix,jx]=float(data[i][j+1])
    return FF

def popKinship_new(D,popnames,outgroup=None,fprefix=None,keep_outgroup=False,hzy=None):
    '''
    Compute kinship matrix between populations

    Parameters:
    ---------------

    D : Reynolds distances between populations
    popnames : names of populations
    outgroup : population to use as outgroup for rooting (if None, use hzy (should give good results) 
               or do midpoint rooting (not implemented))
    fprefix : if set, save the temporary files with names built from fprefix
    keep_outgroup : do not drop the population used for rooting the tree
    hzy : Optimize root finding on observed heterozygosities vector hzy

    Returns:
    -----------

    Fij : Kinship matrix between populations
    fprefix_{ }.txt : if fprefix was set

    See Also:
    ------------

    popgen.reynolds to compute reynolds distances between populations
    popgen.heterozygosity to compute heterozygosities of populations
    '''
    if outgroup is not None:
        try:
            assert outgroup in popnames
        except AssertionError:
            print "Outgroup not found in file"
            print popnames
            sys.exit(1)
    if not fprefix:
        buf,rey_file_name=tempfile.mkstemp(prefix='reynolds')
        buf,fij_file_name=tempfile.mkstemp(prefix='fij',suffix='.txt')
        buf,tree_file_name=tempfile.mkstemp(prefix='tree',suffix='.txt')
    else:
        rey_file_name=fprefix+'_reynolds.txt'
        fij_file_name=fprefix+'_fij.txt'
        tree_file_name=fprefix+'_tree.txt'
    ## write reynolds to file
    reynolds_out=open(rey_file_name,'w')
    for i in range(D.shape[0]):
        tw=[popnames[i]]
        for j in range(D.shape[1]):
            tw.append(str(D[i,j]))
        print >>reynolds_out,' '.join(tw)
    reynolds_out.close()
    ## Build the tree
    dist_mat=2*D
    my_nj=nj.NJ(dist_mat,popnames)
    my_nj.fit()
    if outgroup:
        root_edge=[e for e in my_nj.edges if e.n1.label==outgroup or e.n2.label==outgroup]
        tree=nj.Rooted_Tree()
        tree.build_from_edges(my_nj.edges,root_edge[0])
        if keep_outgroup:
            if hzy is not None:
                ## we can optimize root placement
                ## note that ow the root is in the center of the edge
                hzy_dict={}
                for i,nom in enumerate(popnames):
                    hzy_dict[nom]=hzy[i]
                tree.optim_root(hzy_dict)
        else:
            ## we drop the outgroup
            if tree.root.child[0].label==outgroup:
                tree.root=tree.root.child[1]
            else:
                tree.root=tree.root.child[0]
            tree.root.parent=None
    else:
        ## No outgroup information
        if hzy is not None:
            hzy_dict={}
            for i,nom in enumerate(popnames):
                hzy_dict[nom]=hzy[i]
            ## we can optimize using heterozygosities
            resid=np.inf ## optim criterion
            root_edge_best=None ## where we will place the root
            for candidate in my_nj.edges:
                edges=[e for e in my_nj.edges]
                tree=nj.Rooted_Tree()
                tree.build_from_edges(edges,candidate)
                tree_fit=tree.optim_root(hzy_dict)
                print 'Edge',candidate,' : ',tree_fit.fun
                if tree_fit.x[1]>0 and tree_fit.fun<resid:
                    resid=tree_fit.fun
                    root_edge_best=candidate
                    root_pos_best=tree_fit.x[1]
            tree=nj.Rooted_Tree()
            tree.build_from_edges(my_nj.edges,root_edge_best)
            tree.shift_root(root_pos_best)
        else:
            ## no hzy information, should revert to midpoint but this is
            ## usually not a good idea ...
            raise NotImplementedError
    ## Now we have a tree, get Kinship
    with open(tree_file_name,'w') as f:
        print >>f,tree.newick()
    kin,poplabels=tree.kinship()
    with open(fij_file_name,'w') as f:
        for i,nom in enumerate(poplabels):
            print >>f,nom,' '.join([str(x) for x in kin[i,]])
    ### we have to be careful because the pop order can change within Fij
    pname_temp=popnames[:]
    if outgroup and not keep_outgroup:
        pname_temp.remove(outgroup)
    ordre=[pname_temp.index(x) for x in poplabels]
    FF=kin[ordre,:][:,ordre]
  
    ## clean up temp files if needed
    if not fprefix:
        os.remove(rey_file_name)
        os.remove(fij_file_name)
        os.remove(tree_file_name)
    return FF

def popKinship(D,popnames,outgroup=None,fprefix=None,keep_outgroup=False):
    '''
    Compute kinship matrix between populations

    Parameters:
    ---------------

    D : Reynolds distances between populations
    popnames : names of populations
    outgroup : population to use as outgroup for rooting (if None, do midpoint rooting)
    fprefix : if set, save the temporary files with names built from fprefix
    
    Returns:
    -----------

    Fij : Kinship matrix between populations
    fprefix_{ }.txt : if fprefix was set

    See Also:
    ------------

    popgen.reynolds to compute reynolds distances between populations
    '''
    if outgroup is not None:
        try:
            assert outgroup in popnames
        except AssertionError:
            print "Outgroup not found in file"
            print popnames
            sys.exit(1)
    if not fprefix:
        buf,rey_file_name=tempfile.mkstemp(prefix='reynolds')
        buf,R_file_name=tempfile.mkstemp(prefix='compFij',suffix='.R')
        buf,fij_file_name=tempfile.mkstemp(prefix='fij',suffix='.txt')
        buf,tree_file_name=tempfile.mkstemp(prefix='tree',suffix='.txt')
    else:
        rey_file_name=fprefix+'_reynolds.txt'
        R_file_name=fprefix+'_kinship.R'
        fij_file_name=fprefix+'_fij.txt'
        tree_file_name=fprefix+'_tree.txt'
    ## write reynolds to file
    reynolds_out=open(rey_file_name,'w')
    for i in range(D.shape[0]):
        tw=[popnames[i]]
        for j in range(D.shape[1]):
            tw.append(str(D[i,j]))
        print >>reynolds_out,' '.join(tw)
    reynolds_out.close()
    ## create R script
    R_script=open(R_file_name,'w')
    print >>R_script,Fij_Rscript
    print >>R_script,'D=read.table("'+rey_file_name+'",row.names=1)'
    if outgroup is not None:
        if keep_outgroup:
            print >>R_script,'F=kinship(D,outgroup="'+outgroup+'",keep_outgroup=TRUE)'
        else:
            print >>R_script,'F=kinship(D,outgroup="'+outgroup+'")'
    else:
        print >>R_script,'F=kinship(D)'
    print >>R_script,'write.table(F$matrix,"'+fij_file_name+'",quote=FALSE,col.names=FALSE)'
    if fprefix is not None:
        print >>R_script,'pdf(file="'+fprefix+'_popTree.pdf",w=10,h=10)'
        print >>R_script,'plot(F$tree)'
        print >>R_script,'dev.off()'
        print >>R_script,'write.tree(F$tree,file="'+tree_file_name+'")'
    R_script.close()
    os.system('Rscript --no-save --no-restore '+R_file_name+' 2>&1 1>/dev/null')
    try:
        data=[x.split() for x in open(fij_file_name).readlines()]
    except:
        print '''
Something went wrong when estimating kinship matrix
You need to have R and the 'ape' and 'phangorn' packages installed on your system.
http://cran.r-project.org/web/packages/ape/
http://cran.r-project.org/web/packages/phangorn/
'''
        raise
    ### we have to be careful because the pop order can change within Fij
    pname_temp=popnames[:]
    if outgroup and not keep_outgroup:
        pname_temp.remove(outgroup)
    ordre=[pname_temp.index(x[0]) for x in data]
    FF=np.zeros((len(pname_temp),len(pname_temp)),dtype=float)
    for i in range(len(pname_temp)):
        ix=ordre[i]
        for j in range(len(pname_temp)):
            jx=ordre[j]
            FF[ix,jx]=float(data[i][j+1])
    ## clean up temp files if needed
    if not fprefix:
        os.remove(rey_file_name)
        os.remove(R_file_name)
        os.remove(fij_file_name)
        os.remove(tree_file_name)
    return FF

class FLK_result():
    def __init__(self):
        self.p0=np.nan
        self.val=0
        self.pval=1
        
class FLK_test():
    '''
    FLK test implementation
    '''
    def __init__(self,kinship,diallelic=True):
        '''
        Creates an FLK test

        Parameters:
        ---------------

        kinship : population kinship matrix
        diallelic : if True, the locus is diallelic (e.g. a SNP) o.w. multiallelic
        
        '''
        self.F=kinship
        assert self.F.shape[0]==self.F.shape[1]
        self.dimension=self.F.shape[0]
        self.diallelic=diallelic
        self.invF=np.linalg.inv(self.F)
        self.un = np.ones(self.dimension)
        self.w = np.dot(self.invF,self.un.T)/np.dot(self.un,np.dot(self.invF,self.un.T))
        self.D,self.Q = np.linalg.eigh(self.F)
        if diallelic:
            self.eigen_contrib=self.eigen_contrib_diallelic
        else:
            self.eigen_contrib=self.eigen_contrib_multi
        
    def eigen_contrib_diallelic(self,p,diallelic=True):
        '''
        Computes orthogonalized contributions of each component to the FLK
        statistic, using the spectral decomposition of the variance
        covariance matrix V(p) V=cste*F=cste*Q*D*Q.T where Q is the matrix
        of eigen vectors and D the corresponding diagonal matrix with
        eigen values on the diagonal.

        Parameters:
        ---------------

        p : vector of population allele frequencies at a locus

        Return Value:
        -----------------
        a tuple of :
            - p0hat : estimated allele frequency in the ancestral population
            - Z : contributions of each PC to the FLK statistic
            - cste : multiplying constant to the test
        the value of FLK is = C*sum(L^2)
        '''
        p0hat = np.dot(self.w.T,p)
        cste = ((1 - (1/(np.dot(self.un,np.dot(self.invF,self.un.T)))))/(p0hat*(1-p0hat)))
        Z=np.dot((p-(p0hat*self.un)),self.Q)
        Z=np.dot(Z,np.diag(np.sqrt(1/self.D)))
        Z *= np.sqrt(cste)
        # tmp=np.dot(p-(p0hat*self.un),self.invF*cste)
        # tmp=np.dot(tmp,(p-(p0hat*self.un).T))
        # print tmp,np.sum([z**2 for z in Z])
        # R = np.dot(self.Q.T,(p-(p0hat*self.un).T))
        # R = np.dot(np.diag(np.sqrt(1/self.D)),R)
        # R *= np.sqrt(cste)
        # print R.shape
        # Z=[R[i] for i in range(R.shape[0])]
        return p0hat,Z,cste

    def eigen_contrib_multi(self,p):
        '''
        Multiallelic version of eigen_contrib_diallelic
        
        Parameters:
        ---------------
        
        p : vector of population allele frequencies at a locus
        
        Return Value:
        -----------------
        a tuple of :
        - p0hat : estimated allele frequency in the ancestral population
        - Z : contributions of each PC to the FLK statistic
        - cste : multiplying constant to the test
        the value of FLK is = C*sum(L^2)
        
        See Also:
        ------------
        
        eigen_contrib_diallelic 
        '''
        p0hat = np.dot(self.w.T,p)
        Z=np.dot((p-(p0hat*self.un)),self.Q)
        Z=np.dot(Z,np.diag(np.sqrt(1/self.D)))
        return p0hat,Z,1.0
    
    def eval_flk(self,p):
        '''
        Computes the FLK statistic

        Parameters:
        ---------------

        p : vector of population allele frequencies

        Return Value:
        -----------------
        a list of :
           -- p0
           -- FLK value
           -- corresponding p-value
           -- each contribution of the PC in turn
        '''
        assert p.shape[0]==self.dimension
        p0,Z,C=self.eigen_contrib(p)
        val=sum([x**2 for x in Z])
        pval=1-chi2.cdf(val,self.dimension-1)
        return [p0,val,pval]+Z.tolist()
      
