''' Implementation of the Neighbour Joining Algorithm '''
from __future__ import print_function
import numpy as np
from scipy.optimize import minimize as optim

class Edge():
    def __init__(self,node1,node2,length):
        self.n1=node1
        self.n2=node2
        self.l=length
    def __str__(self):
        return "%s -- [ %f ] -- %s"%(self.n1.label,self.l,self.n2.label)
    
class Node():
    def __init__(self,label=''):
        self.label=label
        self.child=[] ## list of children nodes
        self.f=[]     ## distance to children
        self.parent=None ## parent node
    def newickize(self,parent_distance):
        ''' Create a newick representation of the node and its children'''
        if self.parent is None:
            tw=self.label
        else:
            tw=self.label+':'+str(parent_distance)
        if len(self.child)==0:
            return (tw)
        else:
            subtrees= (sub.newickize(self.f[i]) for i,sub in enumerate(self.child))
            return '(%s)%s' % (','.join(subtrees),tw)       
    def grow_from_edges(self,cur_edges):
        ''' Construct recursively a subtree from edges information '''
        ## we build a new tree: reset values
        self.child=[]
        self.f=[]
        mye=[e for e in cur_edges if e.n1 is self or e.n2 is self]
        if len(mye)==0:
            return
        for e in mye:
            assert e.l>0
            if e.n1 is self:
                self.child.append(e.n2)
                self.f.append(e.l)
                e.n2.parent=self
            elif e.n2 is self:
                self.child.append(e.n1)
                self.f.append(e.l)
                e.n1.parent=self
            cur_edges.remove(e)
        for c in self.child:
            c.grow_from_edges(cur_edges)
        return
    def get_fi(self,curval,res):
        ''' Computes fi : distrance from a tip to the root'''
        if len(self.child)==0:
            ## tip
            res[self.label]=curval
        else:
            for i,c in enumerate(self.child):
                c.get_fi(curval+self.f[i],res)
        return
    def distance_to_root(self,curval=0):
        res=curval
        if self.parent is not None:
            res+=self.parent.f[self.parent.child.index(self)]
            res+=self.parent.distance_to_root(curval)
            return res
        else:
            return 0
    def get_leaves(self):
        ''' Returns a list of leaves (node instances) under self '''
        if len(self.child)==0:
            res=[self]
        else:
            res=[]
            for c in self.child:
                res+=c.get_leaves()
        return res
    def mrca(self,other):
        ''' Identifies most recent common ancestor with other node'''
        try:
            assert len(self.child)==0 and len(other.child)==0
        except AssertionError:
            raise AssertionError('MRCA only applies to leaves')
        anc=self.parent
        while other not in anc.get_leaves():
            anc=anc.parent
            assert anc is not None
        return anc
        
###Optimization function to root the tree 
def qfunc(h0,fi,hi):
    return sum( (((1-fi)*h0-hi)/hi)**2)

def qfunc_2pars(pars,tree,hzy):
    h0=pars[0]
    posroot=pars[1]
    fi_tips=tree.calc_fi(posroot)
    fi_opt=np.zeros(len(fi_tips),dtype=float)
    hi_opt=np.zeros(len(fi_tips),dtype=float)
    for i,pop in enumerate(fi_tips.keys()):
        fi_opt[i]=fi_tips[pop]
        hi_opt[i]=hzy[pop]
    return(qfunc(h0,fi_opt,hi_opt))

class Rooted_Tree():
    def __init__(self):
        self.root=None
    def build_from_edges(self,edges,root_edge):
        ''' Builds a tree from a list of edges, positioning the root in the 
        center of root_edge'''
        current_edges=[e for e in edges if e is not root_edge]
        self.root=Node(label=root_edge.n1.label+'-'+root_edge.n2.label)
        self.root.child=[root_edge.n1,root_edge.n2]
        self.root.f=[0.5*root_edge.l,0.5*root_edge.l]
        for c in self.root.child:
            c.parent=self.root
            c.grow_from_edges(current_edges)
    def newick(self):
        if self.root:
            return self.root.newickize('')+';'
    def shift_root(self,posroot):
        root_edge_len=np.sum(self.root.f)
        self.root.f[0]=posroot*root_edge_len
        self.root.f[1]=(1-posroot)*root_edge_len
    def calc_fi(self,posroot=None):
        '''
        Get a dictionary of distances from tips to root.
        optional parameter posroot (0<pos<1) adjusts the position 
        of the root between its two children
        '''
        if posroot:
            self.shift_root(posroot)
        res={}
        self.root.get_fi(0,res)
        return res
    def optim_h0(self,hzy):
        fi_tips=self.calc_fi()
        fi_opt=np.zeros(len(fi_tips),dtype=float)
        hi_opt=np.zeros(len(fi_tips),dtype=float)
        for i,pop in enumerate(fi_tips.keys()):
            fi_opt[i]=fi_tips[pop]
            hi_opt[i]=hzy[pop]
        val=optim(qfunc,np.mean(hi_opt),args=(fi_opt,hi_opt))
        return val
    def optim_root(self,hzy):
        val=optim(qfunc_2pars,[0.25,0.5],args=(self,hzy))
        return val
    def kinship(self):
        ''' Returns the Kinship matrix of the leaves '''
        leaves=self.root.get_leaves()
        n=len(leaves)
        kinship=np.zeros((n,n),dtype=float)
        for i,l1 in enumerate(leaves):
            kinship[i,i]=l1.distance_to_root()
            for j in range(i+1,n):
                anc=l1.mrca(leaves[j])
                kinship[i,j]=anc.distance_to_root()
                kinship[j,i]=kinship[i,j]
        return kinship,[l.label for l in leaves]
    
class NJ():
    ''' Class Implementing the Neighbour Joining Algorithm
    '''
    def __init__(self,Distance,labels):
        ''' 
        Create an instance from a Distance Matrix n x n and a 
        list of corresponding population albels '''
        assert len(labels)==Distance.shape[0] and Distance.shape[0]==Distance.shape[1]
        self.D=Distance
        self.tips=[Node(x) for x in labels]
        self.Dnodes=[x for x in self.tips]
        self.edges=[]
        
    def _QfromD(self,D):
        ''' Calculates Q-matrix from D.
        Returns Q and (imax,jmax) the indices of smallest Q'''
        n=D.shape[0]
        rowsum=np.sum(D,axis=0)
        Q=(n-2)*D
        for i in range(n):
            Q[:,i]-=rowsum[i]
            Q[i,:]-=rowsum[i]
            Q[i,i]=0
        ## make sure the diagonal does not contain smallest value
        Q+=np.diag(n*[np.max(Q)+1])
        return Q,np.unravel_index(np.argmin(Q),Q.shape)

    def _get_new_D(self,D,pair):
        '''Calculates a new distance matrix after joining the pair (i,j) into new node u 
           and creating the 2 new corresponding edges.
           return newD
        '''
        n=D.shape[0]
        delta=np.zeros(2)
        imin,jmin=pair
        node1=self.Dnodes[imin]
        node2=self.Dnodes[jmin]
        node_intern=Node(label=node1.label+'-'+node2.label)
        delta[0]=0.5*D[imin,jmin]+(0.5/(n-2))*(np.sum(D[imin,])-np.sum(D[jmin,]))
        delta[1]=D[imin,jmin]-delta[0]
        ## create new edges
        new_edge=Edge(node1,node_intern,delta[0] >0 and delta[0] or 0)
        self.edges.append(new_edge)
        new_edge=Edge(node2,node_intern,delta[1]>0 and delta[1] or 1)
        self.edges.append(new_edge)

        newD=np.zeros((n-1,n-1))
        ind_left=range(n)
        ind_left.remove(imin)
        ind_left.remove(jmin)
        self.Dnodes.remove(node1)
        self.Dnodes.remove(node2)
        self.Dnodes=[node_intern]+self.Dnodes
        newD[1:,1:]=D[ind_left,:][:,ind_left]
        for newi,i in enumerate(ind_left):
            newD[0,newi+1]=0.5*(D[imin,i]+D[jmin,i]-D[imin,jmin])
            newD[newi+1,0]=newD[0,newi+1]
        return newD
    
    def fit(self):
        '''
        Calculates the neighbor joining tree
        Sets up internal nodes and edges between nodes.
        '''
        curD=np.copy(self.D)
        while curD.shape[0]>2:
            # print('+++')        
            Q,pair=self._QfromD(curD)
            # print(Q,pair,sep='\n')
            newD=self._get_new_D(curD,pair)
            # print(self.edges[-1])
            # print(self.edges[-2])
            # print(newD,sep='\n')
            curD=newD
        ## terminate
        new_edge=Edge(self.Dnodes[0],self.Dnodes[1],curD[0,1])
        self.edges.append(new_edge)
        # print(self.edges[-1])

       
    
def test():
    ## Test distance matrix from wikipedia :)
    print("#### Wikipedia example ####\nsee https://en.wikipedia.org/wiki/Neighbor_joining")
    Dtest=np.array([[0,5,9,9,8],
                    [5,0,10,10,9],
                    [9,10,0,8,7],
                    [9,10,8,0,3],
                    [8,9,7,3,0]])
    labels=['a','b','c','d','e']
    my_nj=NJ(Dtest,labels)
    my_nj.fit()
    print("Edges from NJ:")
    print(*my_nj.edges,sep='\n')
    print("#### Ovis  example ####")
    Dsheep=2*np.array([[0.0, 0.052524410884, 0.479593743725, 0.15280193534],
                     [0.052524410884 , 0.0 , 0.489188722484 ,0.162238830969],
                     [0.479593743725 , 0.489188722484 , 0.0 , 0.423091249114],
                     [0.15280193534 , 0.162238830969 , 0.423091249114 ,0.0]])
    labels=['MOOA','IROA','Capra','IROO']
    my_nj=NJ(Dsheep,labels)
    my_nj.fit()
    print("Edges from NJ:")
    print(*my_nj.edges,sep='\n')
    ## remove outgroup
    outgroup='Capra'
    ## find edge from root
    root_edge=[e for e in my_nj.edges if e.n1.label==outgroup or e.n2.label==outgroup][0]
    tree=Rooted_Tree()
    tree.build_from_edges(my_nj.edges,root_edge)
    if tree.root.child[0].label==outgroup:
        tree.root=tree.root.child[1]
    else:
        tree.root=tree.root.child[0]
    tree.root.parent=None
    print("Outgroup Tree:",tree.newick())
    Dsheep_nocapra=Dsheep[(0,1,3),:][:,(0,1,3)]
    labels.remove('Capra')
    my_nj=NJ(Dsheep_nocapra,labels)
    my_nj.fit()
    print(*my_nj.edges,sep='\n')
    hzy={'MOOA':0.2531,'IROA':0.2479,'IROO':0.2874}
    resid=np.inf
    root_edge_best=None
    with open('sheep.tree','w') as ftree:
        for root_e in my_nj.edges:
            the_edges=[e for e in my_nj.edges]
            tree=Rooted_Tree()
            tree.build_from_edges(the_edges,root_e)
            print(resid)
            # tree_fit=tree.optim_h0(hzy)
            tree_fit=tree.optim_root(hzy)
            print(tree_fit.x,tree_fit.fun)
            if tree_fit.x[1]>0 and tree_fit.fun<resid:
                resid=tree_fit.fun
                root_edge_best=root_e
                root_pos_best=tree_fit.x[1]
    print("Found best tree with %f accuracy"%resid)
    tree=Rooted_Tree()
    tree.build_from_edges(my_nj.edges,root_edge_best)
    tree.shift_root(root_pos_best)
    kin,popnames=tree.kinship()
    print(kin,*popnames)
    
if __name__=='__main__':
    test()
