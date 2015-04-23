import scipy as sp
import scipy.io
import scipy.optimize as opt
import numpy as np

def MNE_test():
    example = sp.io.loadmat('example_MNE.mat')
    pstart = np.reshape(example['pstart'], (-1))
    resp = example['spikes']
    stim = example['stimuli']

    order = 1;

    Nsamples, Ndim = stim.shape
    psp = np.mean(resp) #spike probability
    avg = stim.T.dot(resp) / Nsamples
    avgs = np.vstack((psp,avg))
    if order > 1 :
        avgsqrd = stim.T.dot(np.tile(resp, (1,Ndim)) * stim) / Nsamples
        avgsqrd = np.reshape(avgsqrd,(Ndim**2,1))
        avgs = np.vstack((avgs,avgsqrd))

    ll_pstart = log_loss(pstart, stim, resp, order)
    dll_pstart = d_log_loss(pstart, stim, avgs, order)

    print 'difference in evaluated log loss'
    print ll_pstart - example['ll_pstart']
    print 'difference in evaluated d log loss'
    print np.max(dll_pstart - example['dll_pstart'])

def log_loss(p, stim, resp, order):
    #get number of samples and dimensionality of stimulus
    Nsamples, Ndim = stim.shape
    
    #unpack p: (var names match names in Fitzgerald paper)
    a = p[0]
    h = p[1:Ndim+1].T
    
    #case: second order calculation --> need J
    if order > 1:
        #reshape J into Ndim x Ndim matrix:
        J = np.reshape(p[Ndim+1:Ndim+1+Ndim**2], (Ndim,Ndim)).T
    
    if order == 1:
        f1 = 1 + np.exp( a + stim.dot(h))
        f0 = 1 + np.exp(-a - stim.dot(h))
    else:
        f1 = 1 + np.exp( a + stim.dot(h) + (np.sum(stim * (stim.dot(J)),1)))
        f0 = 1 + np.exp(-a - stim.dot(h) - (np.sum(stim * (stim.dot(J)),1)))
    
    F1 = resp * np.log(f1)
    F0 = (1 - resp) * np.log(f0)
    F1[np.isnan(F1)] = 0
    F0[np.isnan(F0)] = 0
    return np.mean(F0 + F1)

def d_log_loss(p,stim,avgs,order):
    #get number of samples and dimensionality of stimulus
    Nsamples, Ndim = stim.shape
    
    #unpack p: (var names match names in Fitzgerald paper)
    a = p[0]
    h = p[1:Ndim+1].T
    
    #case: second order calculation --> need J
    if order > 1:
        J = np.reshape(p[Ndim+1:Ndim+1+Ndim**2], (Ndim,Ndim))
        
    if order == 1:
        pSpike = 1.0 / (1.0 + np.exp(a + stim.dot(h))) #Nsamples x 1
        averages = np.hstack((np.mean(pSpike), stim.T.dot(pSpike) / Nsamples))
    elif order == 2:
        
        pSpike = 1.0 / (1.0 + np.exp(a + stim.dot(h) + (np.sum(stim * (stim.dot(J)),1))))
        averages = np.zeros(1+Ndim+Ndim**2)
        averages[0] = np.mean(pSpike)
        averages[1:Ndim+1] = stim.T.dot(pSpike) / Nsamples #ave number of spikes for each stim dimension
        
        temp = (stim.T.dot(np.tile(np.reshape(pSpike, (Nsamples, 1)), (1,Ndim)) * stim)) / Nsamples  #ave number of spikes for each stim correlation
        temp = np.reshape(temp,[Ndim**2,1])
        averages[Ndim+1:Ndim+1+Ndim**2] = np.reshape(temp, Ndim**2)
        
    return (np.squeeze(avgs) - averages)

if __name__ == "__main__":
    MNE_test()