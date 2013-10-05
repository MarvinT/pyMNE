import numpy as np

'''
Log Loss function: 
    L = -<Log[P(y|x)]>, averaged over the empirical distribution
    L = -1/Nsamples * sum_Stim(NspikesperStim/Nrep * log[P(spike|Stim)] + NsilencesperStim/Nrep * log[P(no spike|Stim)])
Log loss function = noise entropy -> maximizing the noise entropy of a joint ensemble yields
the least biased estimate of the conditional probability

arguments: (all of matrix type)
    p: vector where the following holds:
        P(spike|s)^-1 = 1 + exp(p(1) + p(2:Ndim_1)*s + s'*J*s) (note: matlab notation)
        where J is reshaped form of p(Ndim+2:Ndim+1+Ndim^2) 
        ... assume: p = 1 dimensional matrix (== row vector)
        
    stim: (Nsamples x Ndim) every row is a sample == represented by Ndim vector
    
    resp: (Nsamples x 1) == average neuronal population response to each sample
    
    order (1 or 2 (or just greater)) == order of correlation

'''
def log_loss(p, stim, resp, order):
    #get number of samples and dimensionality of stimulus
    Nsamples = np.size(stim,0)
    Ndim = np.size(stim,1)
    
    #unpack p: (var names match names in Fitzgerald paper)
    a = p[0,0]
    h = p[0,1:Ndim+1]
    
    #case: second order calculation --> need J
    if(order > 1):
        J_squash = p[0,Ndim+1:Ndim+2+Ndim**2]
        #reshape J into Ndim x Ndim matrix:
        J = np.reshape(J_squash,(Ndim,Ndim)) #matrix?
    
    if(order == 1):
        f1 = 1 + np.exp(a + stim * h.T)
        f0 = 1 + np.exp(-a - stim * h.T)
    else:
        f1 = 1 + np.exp(np.array(np.ones(Nsamples)*a) + np.array(stim * h.T)[:,0] + (np.sum(np.array(stim)*np.array(stim*J),1)))
        f0 = 1 + np.exp(-np.array(np.ones(Nsamples)*a) - np.array(stim * h.T)[:,0] - (np.sum(np.array(stim)*np.array(stim*J),1)))
        
    F1 = np.array(resp)[:,0] * np.log(np.array(f1))
    F0 = (1 - np.array(resp)[:,0]) * np.log(np.array(f0))
    F1[np.isnan(F1)] = 0
    F0[np.isnan(F0)] = 0
    return np.mean(F0 + F1)
    
def test_log_loss():
    '''
    Nsamples = 3
    Ndim = 2
    '''
    stim = np.matrix([[1,0],[0,1],[1,1]])
    resp = np.matrix([[1],[0],[1]])
    p_order1 = np.matrix([[.5],[.5],[1]]).T
    print log_loss(p_order1,stim,resp,1) #matches matlab code = 1.2139
    p_order2 = np.matrix([[.5],[.5],[.5],[.1],[.2],[.3],[.4]]).T
    print log_loss(p_order2,stim,resp,2) #matches matlab code = 1.3955
    
#test_log_loss()



'''
Gradient of the log loss function:
    moves p-vector towards values that will generate the STA and STC of the
    population response
    
    ...input values are in matrix form
'''
def d_log_loss(p,stim,avgs,order):
    #get number of samples and dimensionality of stimulus
    Nsamples = np.size(stim,0)
    Ndim = np.size(stim,1)
    
    #unpack p: (var names match names in Fitzgerald paper)
    a = p[0,0]
    h = p[0,1:Ndim+1]
    
    #case: second order calculation --> need J
    if(order > 1):
        J_squash = p[0,Ndim+1:Ndim+2+Ndim**2]
        #reshape J into Ndim x Ndim matrix:
        J = np.reshape(J_squash,(Ndim,Ndim))
        
    if(order == 1):
        pSpike = 1/(1 + np.exp(a + stim * h.T)) #Nsamples x 1
        averages = np.zeros(Ndim+1)
        averages[0] = np.mean(pSpike)
        averages[1:] = np.array(stim.T*pSpike)[:,0]/Nsamples #Nsamples x 1
    else: #assume oreder = 2
        pSpike = 1 / (1 + np.exp(np.array(np.ones(Nsamples)*a) + np.array(stim * h.T)[:,0] + (np.sum(np.array(stim)*np.array(stim*J),1))))
        averages = np.zeros(np.size(p))
        averages[0] = np.mean(pSpike)
        averages[1:Ndim+1] = np.array(stim.T*np.matrix(pSpike).T)[:,0]/Nsamples
        temp = (stim.T * (np.array(np.tile(pSpike,(Ndim,1))).T * np.array(stim)))/Nsamples
        temp = np.reshape(temp,[Ndim**2,1])
        averages[Ndim+1:Ndim+1+Ndim**2] = np.array(temp)[:,0]
        
    return np.array(avgs)[:,0] - averages
    

def test_d_log_loss():
    '''
    Nsamples = 3
    Ndim = 2
    '''
    stim = np.matrix([[1,0],[0,1],[1,1]])
    avgs1 = np.matrix([[.1],[.2],[.3]])
    p_order1 = np.matrix([[.5],[.5],[1]]).T
    print d_log_loss(p_order1,stim,avgs1,1) #matlab check: -0.0902    0.0706    0.1995
    p_order2 = np.matrix([[.5],[.5],[.5],[.1],[.2],[.3],[.4]]).T
    avgs2 = np.matrix([[.1],[.2],[.3],[.4],[.5],[.6],[.7]])
    print d_log_loss(p_order2,stim,avgs2,2) 
    # matches matlab code:  -0.0745    0.0915    0.2088    0.2915    0.4747    0.5747    0.6088
    
test_d_log_loss()

    
    
            
    
    
    
    
    
    
    
    