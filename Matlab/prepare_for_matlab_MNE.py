import numpy as np 
from ephys import core, events, spiketrains
from scipy.io import savemat

def create_subwindows(segment, subwin_len):
    ''' Create list of subwindows for cell group identification 

    Parameters
    ------
    segment : list
        Beginning and end of the segment to subdivide into windows
    subwin_len : int
        number of samples to include in a subwindows
    n_subwin_starts : int
        number of shifts of the subwindows

    Returns
    ------
    subwindows : list 
        list of subwindows
    '''
    

    starts = [np.floor(segment[0])]
    subwindows = []
    for start in starts:
        subwin_front = np.arange(start, segment[1], subwin_len)
        for front in subwin_front:
            subwin_end = front + subwin_len
            subwindows.append([front, subwin_end])

    return subwindows

def prepare_for_matlab_MNE(block_path, rec, stimulus, unit, winsize=21.0):
	'''
	Prepares input files for matlab MNE computations from ephys-analysis form 

	Parameters
	------
	block_path : str 
		Path to folder containing post-manually sorted data 
	stimulus : str 
		Name of stimulus you are computing MNE for 
	unit : int 
		Unit ID of unit you wish to compute MNE for
	winsize : float 
		size in milliseconds of the time windows.  Shouldn't need to change 
	'''

	spikes = core.load_spikes(block_path)
	fs = core.load_fs(block_path)
	trials = events.load_trials(block_path)

	# Get all spikes from the unit of interest
	spikes = spikes[spikes['cluster']==unit]

	# Get all trials from the stimulus of interest 
	trials = trials[trials['stimulus']==stimulus]

	stimulus_duration_samples = trials.iloc[0]['stimulus_end'] - trials.iloc[0]['time_samples']
	window_size_samples = np.round(winsize/1000. * fs)
	maxwin = np.floor_divide(stimulus_duration_samples, window_size_samples)
	ntrials = trials.shape[0]

	response_data = np.zeros((ntrials, maxwin-1))

	for trialnum in range(ntrials):
		trial_start = trials.iloc[trialnum]['time_samples']
		trial_end = trials.iloc[trialnum]['stimulus_end']
		# Get spikes in trial
		spike_mask = ((spikes['time_samples']>trial_start) & (spikes['time_samples']<=trial_end))
		spike_times_in_trial = spikes['time_samples'][spike_mask].values - trial_start
		spikes_binned = np.floor_divide(spike_times_in_trial, window_size_samples)
		for win in range(maxwin-1):
			response_data[trialnum, win] = np.sum(1.0*(spikes_binned==win))

	response_matfilename = os.path.join(block_path, 'MNE_response_mat-stim_{}-unit_{}'.format(stimulus, unit))		
	savemat(response_matfilename, {'response':response_data})

def prepare_stimulus_for_MNE(stimulus_file):

	#load wavefile
	(fs, stim_data) = wavefile.read(stimulus_file)

	#downsample
	q = np.round(fs/24000.)
	stim_deci = decimate(stim_data, q, zero_phase=True)
	



