import scipy.optimize as opt
import numpy as np

def log_loss(p, stim, resp, order):
    #get number of samples and dimensionality of stimulus
    Nsamples, Ndim = stim.shape
    resp = np.reshape(resp, (-1))
    
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

def constrained_averages(stim, resp, order):
    Nsamples, Ndim = stim.shape
    psp = np.mean(resp) #spike probability
    avg = stim.T.dot(resp) / Nsamples
    avgs = np.vstack((psp,avg))
    if order > 1 :
        avgsqrd = stim.T.dot(np.tile(resp, (1,Ndim)) * stim) / Nsamples
        avgsqrd = np.reshape(avgsqrd,(Ndim**2,1))
        avgs = np.vstack((avgs,avgsqrd))
    return avgs

def rand_pstart(avgs, order, Ndim):
    pstart = np.log(1.0 / avgs[0] - 1.0)
    pstart = np.hstack((pstart,(.001*(2*np.random.rand(Ndim)-1))))
    if order > 1:
        temp = .0005 * (2 * np.random.rand(Ndim,Ndim) - 1) # for symmetry
        pstart = np.hstack((pstart, np.reshape(temp + temp.T, Ndim**2)))
    return pstart
    
class IterCounter(object):
    def __init__(self):
        self.n_iters = 1
        print '{0:5s}'.format('Iters')
    def callback(self, xk):
        print '{0:5d}'.format(self.n_iters)
        self.n_iters += 1
    
class OverfitException(Exception):
    def __init__(self, p):
        self.p = p
        
class OverfitStopper(object):
    def __init__(self, test_stim, test_resp, order):
        self.test_stim = test_stim
        self.test_resp = test_resp
        self.order = order
        self.n_iters = 1
        self.best_ll = np.inf
        self.best_p = 0
        self.test_tally = 0
        print '{0:5s}   {1:5s}   {2:9s}'.format('Iters', 'tally', 'll(test)')
        
    def callback(self, pk):
        ll_test_k = log_loss(pk, self.test_stim, self.test_resp, self.order)
        print '{0:5d}   {1:5d}   {2: 3.6f}'.format(self.n_iters, self.test_tally, ll_test_k)
        if self.n_iters <= 2 or ll_test_k < self.best_ll:
            self.best_ll = ll_test_k
            self.best_p = pk
            self.test_tally = 0
        else:
            self.test_tally += 1
        
        if self.test_tally >= 10:
            print 'minimum of test set found'
            raise OverfitException(self.best_p)
        
        self.n_iters += 1

def MNEfit(stim, resp, order, pstart=None):

    stim = np.array(stim, dtype=float)
    resp = np.array(resp, dtype=float)

    Nsamples, Ndim = stim.shape
    avgs = constrained_averages(stim, resp, order)
    
    if pstart is None: #initialize params:
        pstart = rand_pstart(avgs, order, Ndim)
    
    #redefine functions with fixed vals:
    def logLoss(p):
        return log_loss(p, stim, resp, order)
    def dlogLoss(p):
        return d_log_loss(p, stim, avgs, order)
    
    pfinal = opt.fmin_cg(logLoss, pstart, fprime=dlogLoss, 
                         callback=IterCounter().callback, maxiter=200)
    
    return pfinal

def MNEfit_jackknives(stim, resp, order, pstart=None, jackknives=4, shuffle=True):

    stim = np.array(stim, dtype=float)
    resp = np.array(resp, dtype=float)
    
    Nsamples, Ndim = stim.shape #TODO: rename Nsamples to n_samples
    assert resp.shape[0] == Nsamples
    assert resp.shape[1] == 1
    
    if shuffle:
        shuffled_indxs = range(Nsamples)
        np.random.shuffle(shuffled_indxs)
        stim = stim[shuffled_indxs,:]
        resp = resp[shuffled_indxs,:]
    
    for jackknife in range(jackknives):
        test_stim = stim[jackknife::jackknives,:]
        test_resp = resp[jackknife::jackknives,:]
        train_stim = stim[np.mod(np.arange(Nsamples)-jackknife, jackknives) != 0,:]
        train_resp = resp[np.mod(np.arange(Nsamples)-jackknife, jackknives) != 0,:]
        
        avgs = constrained_averages(train_stim, train_resp, order)

        if pstart is None: #initialize params:
            pstart = rand_pstart(avgs, order, Ndim)
            
        #redefine functions with fixed vals:
        def logLoss(p):
            return log_loss(p, train_stim, train_resp, order)
        def dlogLoss(p):
            return d_log_loss(p, train_stim, avgs, order)

        try:
            pfinal = opt.fmin_cg(logLoss, pstart, fprime=dlogLoss,
                                 callback=OverfitStopper(test_stim, test_resp, order).callback,
                                 maxiter=200)
        except OverfitException as e:
            pfinal = e.p
        
        if jackknife == 0:
            all_pfinals = np.zeros((jackknives, len(pstart)))
        all_pfinals[jackknife, :] = pfinal
    
    return all_pfinals