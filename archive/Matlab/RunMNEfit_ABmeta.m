function RunMNEfit_AB(cellnum,jack, order, respname, stimname, fittype, Nlags, Ndim, savepath)

%clear all; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open stimulus and arrange according to # of time lags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%id=fopen(stimname,'rb');
%[stimulus]=fread(id,'double');
%fclose(id);


% cd('/Users/andreikozlov/Work/bird/Analysis/sortedData/st715/stims');


%cellnum = 750;
Nf = 30; %number of frequency bands in STRF
Nlags = 20; %times the length of a bin gives the time length of STRF
order   = 2;   % order of MNE model to fit
fittype = 0;   % 0 for regular fitting, 1 for random fitting
njack   = 4;   % # jackknives to run (also determines the size of each jackknives)
%Nlags   = 1;   % # time lags to use in spatiotemporal models (=1 is just a
%spatial model)

% load pylearn2FeatMNEs.mat % load stimuli (P vectors produced by spectrogram function)
% load sixsongs.mat
% stimulus = reshape(data,size(data,1),[]);

%load ..\andreifeatures\songs\sixsongsSpikes.mat
load ~/sixsongsSpikesmeta.mat
stimulus = data2;

[Ndim, Nsamples]=size(stimulus);
stimulus=stimulus';
% The stimulus units is changed to be z-score
% (zero mean, normalized by standard deviation),
% but could be put into your own unit system (or 
% maybe your stimulus file is already saved
% in the proper units).
 stimulus = stimulus-repmat(mean(stimulus),[Nsamples,1]);
stim = stimulus./repmat(std(stimulus),[Nsamples,1]);
clear data2 stimulus

% redefine stimulus to include time lags
% if Nlags>1
%     Nsamples = Nsamples - (Nlags-1); %total length of stimulus minus 19 time bins
%     Ndimtotal = Ndim*Nlags; %16x20
%     stim = zeros(Nsamples,Ndimtotal);
%     for i=1:Nlags
%         stim(:,Ndim*(i-1)+1:Ndim*i) = stimulus(i:Nsamples+i-1,:);
%     end
% else
%     stim = stimulus;
% end
% clear stimulus;

Nd = sqrt(Ndim); % # pixels per side of image


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open responses and cut according to # of time lags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% response = double(reshape(spikes(cellnum,:,:),1,[]));
response = double(spikes(cellnum,:));
resp = response';
clear spikes

% Normalize spike counts in each bin by the maximal value for the minimal models (for MID they are not
% normalized
% response = response./max(response);
% resp = response;
% resp = resp(Nlags:length(resp));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% break into training and test sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
njack

pfinal=zeros(njack,1+Nf*Nlags*(1+Nf*Nlags));
pfinalR=zeros(njack,1+Nf*Nlags*(1+Nf*Nlags));
featuresP=[];
featuresPR=[];
featuresN=[];
featuresNR=[];

for jack = 1:njack; %loop over all njacks

Ntest = floor(Nsamples/njack);  % could be changed to do different # of jackknives
teststim = stim(1+(jack-1)*Ntest:jack*Ntest,:);
testresp = resp(1+(jack-1)*Ntest:jack*Ntest);
stim(1+(jack-1)*Ntest:jack*Ntest,:) = [];
resp(1+(jack-1)*Ntest:jack*Ntest) = [];


disp('Starting optimization');
tic
celltype = '';  % ignore this
%MNEfit(stim, resp, teststim, testresp, celltype, cellnum, jack, order, Nd, fittype);

%this is from the MNEfit function:

[Nsamples,Ndim] = size(stim);
psp = mean(mean(resp));   % spike probability
avg = (stim'*resp)/Nsamples;  % Ndim x Nrep
avg = mean(avg,2);  % Ndim x 1
avgs = [psp;avg];
if order>1
    avgsqrd = stim'*(repmat(resp,[1,Ndim]).*stim)/Nsamples;  % Ndim x Ndim
    avgsqrd = reshape(avgsqrd,[Ndim^2,1]);
    avgs = [avgs;avgsqrd];
end
% Initialize parameters
pstart = log(1/avgs(1)-1);
pstart(2:Ndim+1) = .001*(2*rand([1,Ndim])-1);
if order>1
    temp = .0005*(2*rand([Ndim,Ndim])-1);
    pstart(Ndim+2:length(pstart)+Ndim^2) = reshape((temp+temp'),[1,Ndim^2]);
    clear temp;
end

% Run conjugate gradient algorithm
pfinal(jack,:) = frprmn(pstart, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, Nd, fittype);
disp(['Optimization took ' num2str(toc/60) ' minutes']);
tic
pfinalR(jack,:) = frprmn(pstart, @logloss, @dlogloss, stim, resp(randperm(size(resp,1))), teststim, testresp(randperm(size(testresp,1))), order, avgs, Nd, fittype);
disp(['Optimization took ' num2str(toc/60) ' minutes']);


%Plot the results and save the figures

% 
% h=pfinal(2:Nlags*Nf+1);
% 
% 
%  J=pfinal(Nlags*Nf+2:end); % this is the covariance matrix J
%  [V,D] = eig(reshape(J,Nlags*Nf,Nlags*Nf));
% 
% 
% figure
% subplot(3,3,1)
%  imagesc(reshape(h,Nf,Nlags))
%  axis xy
%  
% subplot(3,3,2)
% imagesc(reshape(feats(:,cellnum),Nf,Nlags))
% axis xy
% axis square
% 
% subplot(3,3,3)
% eigenvalues = diag(D);
% [eigenvalues_sorted,index] = sort(eigenvalues);
% plot(eigenvalues_sorted,'o');
% 
% subplot(3,3,4)
% eig_sorted_1 = V(:,index(1));
% imagesc(reshape(eig_sorted_1,Nf,Nlags))
% axis xy
% 
% subplot(3,3,5)
% eig_sorted_2=V(:,index(2));
% imagesc(reshape(eig_sorted_2,Nf,Nlags))
% axis xy
% 
% subplot(3,3,6)
% eig_sorted_3=V(:,index(3));
% imagesc(reshape(eig_sorted_3,Nf,Nlags))
% axis xy
% 
% subplot(3,3,7)
% eig_sorted_end=V(:,index(end));
% imagesc(reshape(eig_sorted_end,Nf,Nlags))
% axis xy
% 
% subplot(3,3,8)
% eig_sorted_end_1=V(:,index(end-1));
% imagesc(reshape(eig_sorted_end_1,Nf,Nlags))
% axis xy
% 
% subplot(3,3,9)
% eig_sorted_end_2=V(:,index(end-2));
% imagesc(reshape(eig_sorted_end_2,Nf,Nlags))
% axis xy
% 
% 
% 
% figure_name = (['st_' num2str(cellnum) 'your figure name' '_Nlags' num2str(Nlags) '_nfft128_Nf16' '_jack_' num2str(jack) '_of_' num2str(njack)]); % the 6 songs
% 
% hgsave(figure_name); 


end
[Jreal,Jrand,hreal,hrand]=testSignificance(pfinal,pfinalR,Nf,Nlags);
save(['metacell_' num2str(cellnum) '_Nlags_' num2str(Nlags) '_Nf_' num2str(Nf) '.mat'],'pfinal','pfinalR','Jreal','Jrand','hreal','hrand');

end
