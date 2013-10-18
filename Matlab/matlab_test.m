clear all
close all

options.njack = 4;
options.fittype = 0;   % 0 for regular fitting, 1 for random fitting
options.redo_fit = false;

data_loaded = false;
display('loading data')
load spikes.dat
load stimuli.dat
load rf.dat
display('data loaded')

Nsamples = length(spikes);
Ndims = size(stimuli, 1);
pfinals = zeros(options.njack, Ndims + 1);
for jack=1:options.njack
    filename = sprintf('ModelCell_1_MNE1_jack_%i_freq_band_step_1.mat', jack);
    if isempty(dir(filename)) || options.redo_fit
        stim = stimuli';
        resp = spikes';
        Ntest = floor(Nsamples/options.njack);
        teststim = stim(1+(jack-1)*Ntest:jack*Ntest,:);
        testresp = resp(1+(jack-1)*Ntest:jack*Ntest);
        stim(1+(jack-1)*Ntest:jack*Ntest,:) = [];
        resp(1+(jack-1)*Ntest:jack*Ntest) = [];
        
        pfinals(jack,:) = MNEfit(stim, resp, teststim, testresp, 'model', 1, jack, 1, 'this argument is unecessary?', options.fittype, 1);
    else
        load(filename)
        pfinals(jack,:) = pfinal;
    end
end

for i = 1:options.njack
    subplot(2,options.njack, i + (i > (options.njack / 2)) * options.njack / 2)
    imshow(reshape(pfinals(i,2:end), sqrt(Ndims), sqrt(Ndims)), [-.25, 1])
end

subplot(122)
imshow(rf, [-.25, 1])