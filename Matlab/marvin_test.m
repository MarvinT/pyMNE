clear all 
close all

load('Y:\example_MNE.mat')

order = 2;
Nd = 0;
fittype = 1;

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



pfinal = frprmn(pstart, @logloss, @dlogloss, stimuli, spikes, stimuli, spikes, order, avgs, Nd, fittype);

a = pfinal(1);
h = pfinal(2:Ndim+1);
if order>1
    J = reshape(pfinal(Ndim+2:Ndim+1+Ndim^2),[Ndim,Ndim])';
end

%%

figure(1)
subplot(121)
imagesc(reshape(v, 16, 16))
title('original')
axis equal square

subplot(122)
imagesc(reshape(h, 16, 16))
title('MNE - h')
axis equal square

colormap(gray) 

%%
if order > 1
    [evecs,evals]=eig(J);
    [EV,inds] = sort((diag(evals)));
    figure(2)
    imagesc(reshape(eigvec(:,inds(1)), 16, 16))
    axis equal square
    colormap(gray) 
end