clear all 
close all

MNEtype = 1;

if MNEtype == 1
    load('Y:\example_MNE.mat')
elseif MNEtype == 2
    load('Y:\Andre_MNE.mat')
    py_pfinal = pfinal;
elseif MNEtype == 3
    load('Y:\GitHub\pyMNE\Matlab\bt_gendata.mat')
    stimuli = stimuli';
    spikes = spikes';
end

order = 1;
Nd = 0;
fittype = 0;

[Nsamples,Ndim] = size(stimuli);
psp = mean(mean(spikes));   % spike probability
avg = (stimuli'*spikes)/Nsamples;  % Ndim x Nrep
avg = mean(avg,2);  % Ndim x 1
avgs = [psp;avg];
if order>1
    avgsqrd = stimuli'*(repmat(spikes,[1,Ndim]).*stimuli)/Nsamples;  % Ndim x Ndim
    avgsqrd = reshape(avgsqrd,[Ndim^2,1]);
    avgs = [avgs;avgsqrd];
end

if order>1 % make new pstart
    % Initialize parameters
    pstart = log(1/avgs(1)-1);
    pstart(2:Ndim+1) = .001*(2*rand([1,Ndim])-1);
    if order>1
        temp = .0005*(2*rand([Ndim,Ndim])-1);
        pstart(Ndim+2:length(pstart)+Ndim^2) = reshape((temp+temp'),[1,Ndim^2]);
        clear temp;
    end
end

%%
pfinal = frprmn(pstart, @logloss, @dlogloss, stimuli, spikes, stimuli, spikes, order, avgs, Nd, fittype);

%%
pfinal_fminunc = fminunc(@(x) logloss(x, stimuli, spikes, order), pstart);

%%
a = pfinal(1);
h = pfinal(2:Ndim+1);
if order>1
    J = reshape(pfinal(Ndim+2:Ndim+1+Ndim^2),[Ndim,Ndim])';
end

%%

figure(1)
if MNEtype == 1
    subplot(121)
    imagesc(reshape(v, 16, 16))
    title('original')
    axis equal square

    subplot(122)
    imagesc(reshape(h, 16, 16))
elseif MNEtype == 2
    imagesc(reshape(h, 11, 20))
elseif MNEtype == 3
    imagesc(reshape(h, 16, 16))
end
title('MNE - h')
axis equal square

colormap(gray) 

%%
if order > 1
    [evecs,evals]=eig(J);
    [EV,inds] = sort((diag(evals)));
    figure(2)
    if MNEtype == 1 || MNEtype == 3
        imagesc(reshape(evecs(:,inds(1)), 16, 16))
    elseif MNEtype == 2
        imagesc(reshape(evecs(:,inds(1)), 11, 20))
    end
    axis equal square
    colormap(gray) 
    figure(3)
    plot(EV)
    
end

%%
if MNEtype == 3
    figure(4)
    k=1; kmin=0; kmax=100; hk=loop_slider_n(k,kmin,kmax,1);
    while true
        if ~ishandle(hk)
            break
        end
        k = get(hk, 'Value');

        imagesc(reshape(evecs(:,inds(1)), 16, 16).*k/100.0 + reshape(h, 16, 16).*(100.0-k)/100.0)
        colormap(gray)
        uiwait;
    end
end

%%
fstart = logloss(pstart, stimuli, spikes, order);
dfstart = dlogloss(pstart, stimuli, avgs, order)';
%ffinal = logloss(pfinal, stimuli, spikes, order)
%dffinal = dlogloss(pfinal, stimuli, avgs, order)'