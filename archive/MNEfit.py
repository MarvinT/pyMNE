import scipy as sp
import scipy.optimize as opt
import logLossFuncs as LLF

def MNEfit(stim,resp,order):
    # in order for dlogloss to work, we need to know -<g(yt(n),xt)>data
    # == calculate the constrained averages over the data set
    Nsamples = sp.size(stim,0)
    Ndim = sp.size(stim,1)
    psp = sp.mean(sp.mean(resp)) #spike probability (first constraint)
    avg = (1.0*stim.T*resp)/(Nsamples*1.0)
    avgs = sp.vstack((psp,avg))
    if(order > 1):
        avgsqrd = (stim.T*1.0)*(sp.array(sp.tile(resp,(1,Ndim)))*sp.array(stim))/(Nsamples*1.0)
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
    #pfinal = opt.fmin_tnc(logLoss,pstart,fprime=dlogLoss)
    # conjugate-gradient:
    pfinal = opt.fmin_cg(logLoss,pstart,fprime=dlogLoss)
    #pfinal = opt.fmin(logLoss,pstart,fprime=dlogLoss)
    return pfinal
    
def testMNEfit():
    #essentially just random numbers == not a real test of func
    stim = sp.matrix([[0,0],[1,0],[0,1]])
    resp = sp.matrix([[0],[1],[1]])
    #stim = sp.matrix([[0.0],[1.0],[1.0],[0.0]])
    #resp = sp.matrix([[1.0],[0.0],[0.0],[1.0]])
    #order = 1
    #print 'First order MNE fit:'
    #p =  MNEfit(stim,resp,order)
    #print p
    order = 2
    print 'Second order MNE fit:'
    p = MNEfit(stim,resp,order)
    print p
    J = sp.matrix(sp.reshape(p[-4:],(2,2)))
    print J
    h = sp.matrix(p[1:3])
    print h
    #print h*sp.matrix([1,1]).T + sp.matrix([1,1]) * J * sp.matrix([1,1]).T 
    # calculate likelihood of spike:
    print 'spike prob:'
    print 1/(1 + sp.exp(p[0] + h*sp.matrix([1,1]).T + sp.matrix([1,1]) * J * sp.matrix([1,1]).T ))
    
    
testMNEfit()
    