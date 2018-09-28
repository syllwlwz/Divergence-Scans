import sys
import multiprocessing
import numpy as np
import fastphaseCythonMT as fphMT
import fastphaseCython as fph

def global_distance(from_vec,to_vec):
    M=np.zeros((from_vec.shape[0],to_vec.shape[0]))
    for fk in range(from_vec.shape[0]):
        for tk in range(to_vec.shape[0]):
            M[fk,tk]=np.sum((to_vec[tk,]-from_vec[fk,])**2)
    return np.round(M,3)

def greedy_best_match(M):
    n=np.min(M.shape)
    col=range(M.shape[1])
    row=range(M.shape[0])
    rez={}
    cost=0
    permut=-np.ones(n)
    for i in range(n):
        # print row
        my_arg_min=np.argmin(M)
        my_min=np.min(M)
        k1,k2=np.unravel_index(my_arg_min,M.shape)
        cost+=my_min
        rez[(k1,k2)]=my_min
        # print row[k1],col[k2]
        permut[row[k1]]=col[k2]
        del row[k1]
        del col[k2]
        M=np.delete(M,k1,0)
        if n>1:
            M=np.delete(M,k2,1)
    return cost,permut

def fastphase_windowed_optim(K,wsize=20,n_local_em=20,hseq=None,gseq=None,nthread=1,verbose=False):
    '''
    Find a good starting point for the EM algorithm used to estimate a fastphase model
    by fitting independent model on contiguous windows and merging them

    Parameters:
    -----------

    K : number of fastphase clusters
    wsize : size of windows (in number of SNPs)
    n_local_em : number of EM runs to perform on each window
    hseq : list of haplotypes
    gseq : list of genotypes
    nthread : number of CPU to use

    Return Values:
    --------------

    thetas = a starting value for theta parameters in the model
    '''
    if nthread==1:
        myfph=fph
    elif nthread>1:
        myfph=fphMT
    
    L=len(seq[0])
    parset=[]
    thetas=[]
    for iwin,sbeg in enumerate(range(L)[::wsize]):
        send= (sbeg+wsize)> L and L or sbeg+wsize
        my_mod=myfph.fastphase(send-sbeg)
        for i,s in enumerate(seq):
            my_mod.addHaplotype(str(i),s[sbeg:send])
        loglike=-np.inf
        save_par=None
        for e in range(n_local_em):
            if verbose:
                sys.stderr.write("window %2d / %d : em %2d / %d\r"%(iwin+1,np.round(L/wsize),e+1,n_local_em))
                sys.stderr.flush()
            par=my_mod.fit(nClus=K,nthread=nthread)
            if par.loglike>loglike:
                save_par=par
                loglike=par.loglike
        parset.append(save_par)
        ## label switching correction
        my_imp=my_mod.impute([save_par])
        spz=np.zeros((len(seq),K))
        fpz_tmp=np.zeros((len(seq),K))
        for i in range(len(seq)):
            spz[i,]=my_imp[str(i)][1][0][0,]
            fpz_tmp[i,]=my_imp[str(i)][1][0][-1,]
        if iwin==0:
            for ik in range(K):
                thetas.append(parset[0].theta[:,ik])
            kidx=0
        else:
            M=global_distance(fpz.T,spz.T)
            greed_res=greedy_best_match(M)
            kidx=greed_res[1][kidx]
            for ik in range(K):
                thetas[ik]=np.concatenate([thetas[ik],parset[iwin].theta[:,kidx]])
        fpz=fpz_tmp
    return np.array(thetas)

def fastphase_combined_optim(model,par_list,wsize,nthread=1,verbose=False):
    '''
    combine different fit of the EM algorithm into a boosted one

    for now works only on haplotypes
    '''
    L=model.nLoci
    if nthread==1:
        myfph=fph
    elif nthread>1:
        myfph=fphMT
    K=par_list[0].nClus
    thetas=[]
    parset=[]
    nhaps=len(model.haplotypes)
    ngen=len(model.genotypes)
    nobs=ngen+nhaps
    win_list=range(L)[::wsize]
    fdebug=open('fphoptim.debug','a')
    ##print >>fdebug,'K window EM loglike' 
    for iwin,sbeg in enumerate(win_list):
        send= (sbeg+wsize)> L and L or sbeg+wsize
        my_mod=myfph.fastphase(send-sbeg)
        for k,v in model.haplotypes.items():
            my_mod.addHaplotype(k,v[sbeg:send])
        for k,v in model.genotypes.items():
            my_mod.addGenotype(k,v[sbeg:send])
        # find best local EM
        loglike=-np.inf
        save_par=None
        for e,gpar in enumerate(par_list):
            if verbose:
                sys.stderr.write("window %2d / %d : em %2d / %d\r"%(iwin+1,len(win_list),e+1,len(par_list)))
                sys.stderr.flush()
            loc_par=myfph.modParams(send-sbeg,K)
            loc_par.theta=np.copy(gpar.theta[sbeg:send,])
            loc_par.alpha=np.copy(gpar.alpha[sbeg:send,])
            loc_par.rho=np.copy(gpar.rho[sbeg:send,])
            tmp_par=my_mod.fit(nClus=K,params=loc_par,nthread=6,nstep=1)
            print >>fdebug,K,iwin,e,tmp_par.loglike
            if tmp_par.loglike>loglike:
                save_par=loc_par
                loglike=tmp_par.loglike
        parset.append(save_par)
        ## label switching correction
        my_imp=my_mod.impute([save_par],nthread=nthread)
        spz=np.zeros((nobs,K))
        fpz_tmp=np.zeros((nobs,K))
        for i,pair in enumerate(my_mod.haplotypes.items()):
            k,v=pair
            spz[i,]=my_imp[k][1][0][0,]
            fpz_tmp[i,]=my_imp[k][1][0][-1,]
        for i,pair in enumerate(my_mod.genotypes.items()):
            k,v=pair
            probZ=my_imp[k][1][0]
            a=0.5*(np.sum(probZ,axis=1)+np.sum(probZ,axis=2))
            spz[i+nhaps,]=a[0,]
            fpz_tmp[i+nhaps,]=a[-1,]
        if iwin==0:
            for ik in range(K):
                thetas.append(parset[0].theta[:,ik])
            kidx=0
        else:
            M=global_distance(fpz.T,spz.T)
            greed_res=greedy_best_match(M)
            kidx=greed_res[1][kidx]
            for ik in range(K):
                thetas[ik]=np.concatenate([thetas[ik],parset[iwin].theta[:,kidx]])
        fpz=fpz_tmp
    return parset,np.array(thetas).T

def fastphase_windowed_loglike(model,par_list,wsize,verbose=False):
    wlike=[]
    L=model.nLoci
    K=par_list[0].nClus
    nhaps=len(model.haplotypes)
    ngen=len(model.genotypes)
    nobs=ngen+nhaps
    myfph=fph
    win_list=range(L)[::wsize]
    ##print >>fdebug,'K window EM loglike' 
    for iwin,sbeg in enumerate(win_list):
        send= (sbeg+wsize)> L and L or sbeg+wsize
        my_mod=myfph.fastphase(send-sbeg)
        for k,v in model.haplotypes.items():
            my_mod.addHaplotype(k,v[sbeg:send])
        for k,v in model.genotypes.items():
            my_mod.addGenotype(k,v[sbeg:send])
        # find best local EM
        loglike=-np.inf
        save_par=None
        for e,gpar in enumerate(par_list):
            loc_par=myfph.modParams(send-sbeg,K)
            loc_par.theta=gpar.theta[sbeg:send,]
            loc_par.alpha=gpar.alpha[sbeg:send,]
            loc_par.rho=gpar.rho[sbeg:send,]
            tmp_par=my_mod.fit(nClus=K,params=loc_par,nthread=6,nstep=1)
            wlike.append(tmp_par.loglike)
    return wlike
