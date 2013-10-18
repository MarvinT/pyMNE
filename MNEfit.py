import scipy as sp
import scipy.optimize as opt
import logLossFuncs as LLF

def MNEfit(stim,resp,order):
    # in order for dlogloss to work, we need to know -<g(yt(n),xt)>data
    # == calculate the constrained averages over the data set
    Nsamples = sp.size(stim,0)
    Ndim = sp.size(stim,1)
    psp = sp.mean(sp.mean(resp)) #spike probability (first constaint)
    avg = (stim.T*resp)/Nsamples
    avgs = sp.vstack((psp,avg))
    if(order > 1):
        avgsqrd = stim.T*(sp.array(sp.tile(resp,(1,Ndim)))*sp.array(stim))/Nsamples
        avgsqrd = sp.reshape(avgsqrd,(Ndim**2,1))
        avgs = sp.vstack((avgs,avgsqrd))
    #initialize params:
    pstart = sp.log(1/avgs[0,0] - 1)
    pstart = sp.hstack((pstart,(.001*(2*sp.random.rand(Ndim)-1))))
    if(order > 1):
        temp = .0005*(2*sp.random.rand(Ndim,Ndim)-1)
        pstart = sp.hstack((pstart,sp.reshape(temp+temp.T,(1,Ndim**2))[0]))
    
    #redefine functions with fixed vals:
    def logLoss(p):
        return LLF.log_loss(p, stim, resp, order)
    def dlogLoss(p):
        return LLF.d_log_loss(p, stim, avgs, order)
    #run the function:
    pfinal = opt.fmin_tnc(logLoss,pstart,fprime=dlogLoss)
    return pfinal
    
def testMNEfit():
    #essentially just random numbers == not a real test of func
    stim = sp.matrix([[0,0],[0,0],[1,1]])
    resp = sp.matrix([[0],[0],[1]])
    order = 1
    print 'First order MNE fit:'
    print MNEfit(stim,resp,order)
    order = 2
    print 'Second order MNE fit:'
    print MNEfit(stim,resp,order)
    
testMNEfit()
    