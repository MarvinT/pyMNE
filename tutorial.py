import numpy as np
import matplotlib.cm as cm
import pylab

def main():
    # generate white noise stimuli:
    #   consider 2^14 16x16 white noise patches, 
    #   reshaped to form column vector of length 256
    stimuli = np.random.randn(256, 2**14)
    # create gabor rf
    gb = gabor_fn(theta=np.pi/4)
    pylab.imshow(gb, cmap=cm.gray, interpolation="nearest")
    pylab.show()
    np.savetxt('rf.dat', gb)
    # generate spike train
    v = np.resize(gb, (1, np.size(gb)))
    # use 1/(1+exp(v*x+c)) to generate spike train
    p = 1 / (1 + np.exp(v.dot(stimuli)))
    spikes = np.random.rand(np.size(p)) > p
    # save to compare results in matlab
    np.savetxt('stimuli.dat', stimuli)
    np.savetxt('spikes.dat', spikes, fmt='%d')
    # run MNE
    # display recovered rf to generative rf


def gabor_fn(bw=3,gamma=0.5,psi=0,lambd=6,theta=0):
    # bw    = bandwidth, (1)
    # gamma = aspect ratio, (0.5)
    # psi   = phase shift, (0)
    # lambd = wave length, (>=2)
    # theta = angle in rad, [0 pi)
 
    sigma = lambd/np.pi*np.sqrt(np.log(2)/2)*(2**bw+1)/(2**bw-1);
    sigma_x = sigma;
    sigma_y = sigma/gamma;

    sz = 16;

    x = range(np.floor(-sz/2).astype(int), np.floor(sz/2).astype(int))
    y = range(np.floor(sz/2).astype(int), np.floor(-sz/2).astype(int), -1)
    xv, yv = np.meshgrid(x, y)
    # x (right +)
    # y (up +)

    # Rotation 
    x_theta=xv*np.cos(theta)+yv*np.sin(theta);
    y_theta=-xv*np.sin(theta)+yv*np.cos(theta);

    gb=np.exp(-0.5*(x_theta**2/sigma_x**2+y_theta**2/sigma_y**2))*np.cos(2*np.pi/lambd*x_theta+psi);
    return gb

if __name__ == '__main__':
    main()