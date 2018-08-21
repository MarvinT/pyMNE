function RunMNEfit(cellnum,jack, order, respname, stimname, fittype, Nlags, Ndim, savepath)

%clear all; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open stimulus and arrange according to # of time lags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%id=fopen(stimname,'rb');
%[stimulus]=fread(id,'double');
%fclose(id);


% cd('/Users/andreikozlov/Work/bird/Analysis/sortedData/st715/stims');


cellnum = 750;
Nf = 32; %number of frequency bands in STRF
Nlags = 20; %times the length of a bin gives the time length of STRF
order   = 2;   % order of MNE model to fit
fittype = 0;   % 0 for regular fitting, 1 for random fitting
njack   = 4;   % # jackknives to run (also determines the size of each jackknives)
%Nlags   = 1;   % # time lags to use in spatiotemporal models (=1 is just a
%spatial model)

load MNEprepstimspikes.mat % load stimuli (P vectors produced by spectrogram function)
load sixsongs.mat

% stimulus = [P_m_211_s_08_song_60sec_24kHz P_m_211_s_16_song_60sec_24kHz P_m_211_s_17_song_60sec_24kHz P_m_211_s_25_song_60sec_24kHz P_m_462_s2_4_9_24kHz P_m_497_17_s1_24kHz];

stimulus = reshape(data,size(data,1),[]);

% %do this to average time bins
% stimulus_odd_columns = stimulus(:,1:2:end-1);
% %stimulus_odd_columns = stimulus(:,1:2:end);
% stimulus_even_columns = stimulus(:,2:2:end);
% stimulus_averaged = (stimulus_odd_columns+stimulus_even_columns)./2;
% stimulus = stimulus_averaged;
% 
% %and again:
% stimulus_odd_columns = stimulus(:,1:2:end-1);
% %stimulus_odd_columns = stimulus(:,1:2:end);
% stimulus_even_columns = stimulus(:,2:2:end);
% stimulus_averaged = (stimulus_odd_columns+stimulus_even_columns)./2;
% stimulus = stimulus_averaged;
% 
% %and again
% stimulus_odd_columns = stimulus(:,1:2:end-1);
% %stimulus_odd_columns = stimulus(:,1:2:end);
% stimulus_even_columns = stimulus(:,2:2:end);
% stimulus_averaged = (stimulus_odd_columns+stimulus_even_columns)./2;
% stimulus = stimulus_averaged;
% 
% 
% 
% %do this to average adjacent frequency bins (remember that the first one is
% %the  DC component.)
%  stimulus_odd_rows = stimulus(3:2:end,:);
%  stimulus_even_rows = stimulus(2:2:end,:);
%  stimulus = (stimulus_odd_rows+stimulus_even_rows)./2;
% 
%  %and now do this again to average more; note that we begin at 1, not at 3,
%  %because the DC was already removed. Starting with nfft 128, which gives 65
%  %freq. bands, we get 16 frequencies
% 
% stimulus_odd_rows = stimulus(1:2:end,:);
% stimulus_even_rows = stimulus(2:2:end,:);
% stimulus = (stimulus_odd_rows+stimulus_even_rows)./2;
%  
% 
% stimulus = 20*log10(stimulus); %or should it be 10*log10? probably doesn't matter


%And now convert and save as binary for MID analysis
 
%  filename = sprintf('st793_chan_25_26.raw');  
%  fident=fopen(filename,'wb');
%  fwrite(fident,stimulus,'double');
%  fclose(fident);



[Ndim, Nsamples]=size(stimulus);
stimulus=stimulus';
% The stimulus units is changed to be z-score
% (zero mean, normalized by standard deviation),
% but could be put into your own unit system (or 
% maybe your stimulus file is already saved
% in the proper units).
 stimulus = stimulus-repmat(mean(stimulus),[Nsamples,1]);
stimulus = stimulus./repmat(std(stimulus),[Nsamples,1]);


% redefine stimulus to include time lags
if Nlags>1
    Nsamples = Nsamples - (Nlags-1); %total length of stimulus minus 19 time bins
    Ndimtotal = Ndim*Nlags; %16x20
    stim = zeros(Nsamples,Ndimtotal);
    for i=1:Nlags
        stim(:,Ndim*(i-1)+1:Ndim*i) = stimulus(i:Nsamples+i-1,:);
    end
else
    stim = stimulus;
end
%clear stimulus;
Nd = sqrt(Ndim); % # pixels per side of image


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open responses and cut according to # of time lags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fid = fopen(respname,'rb');
%resp = fread(fid,'uint8');
%fclose(fid);



% cd('/Users/andreikozlov/Work/bird/Analysis/sortedData/st750/sortedData/chan_21_22');
% load chan_21_22.mat;



% response = [spike_vector_m_211_s_08_song_60sec; spike_vector_m_211_s_16_song_60sec; spike_vector_m_211_s_17_song_60sec; spike_vector_m_211_s_25_song_60sec; spike_vector_m_462_s2prime; spike_vector_m_497_17_s1];

response = double(reshape(spikest(cellnum,:,:),1,[]));

% response = response'; % this depends on how the response file was
% %prepared:check the rows-columns
% 
% 
% %do this to average time bins pair-wise
% response_odd_columns = response(:,1:2:end-1); 
% %response_odd_columns = response(:,1:2:end);
% response_even_columns = response(:,2:2:end);
% response_averaged = (response_odd_columns+response_even_columns)./2;
% response = response_averaged;
% 
% %and again:
% response_odd_columns = response(:,1:2:end-1);
% %response_odd_columns = response(:,1:2:end);
% response_even_columns = response(:,2:2:end);
% response_averaged = (response_odd_columns+response_even_columns)./2;
% response = response_averaged;
% 
% %and again
% response_odd_columns = response(:,1:2:end-1); 
% %response_odd_columns = response(:,1:2:end); 
% response_even_columns = response(:,2:2:end);
% response_averaged = (response_odd_columns+response_even_columns)./2;
% response = response_averaged;
% 
% 
  response = response';
 
 % Convert to binary and save for MID analysis
 
 
 %filename = sprintf('st793_chan_25_26.isk');
% fident=fopen(filename,'wt');
% % %fident=fopen(filename,'wb');
% fprintf(fident,'%u \n', response)
% %fwrite(fident,response,'uint8');
% fclose(fident);


% Normalize spike counts in each bin by the maximal value for the minimal models (for MID they are not
% normalized

response = response./max(response);


resp = response;
resp = resp(Nlags:length(resp));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% break into training and test sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
pfinal = frprmn(pstart, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, Nd, fittype);


% % Save results
if strcmp(celltype, 'RGC')==1 || strcmp(celltype, 'LGN')==1
    if fittype==0
        save([celltype '_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '.mat'],'pfinal');
    else
        save([celltype '_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '_random.mat'],'pfinal');
    end
else
    if fittype==0
        

   save(['st_' num2str(cellnum) '_m_211_s_08_s16_s17_s25_24kHz' '_Nlags' num2str(Nlags) '_nfft128_Nf16' '_jack_' num2str(jack) '_of_' num2str(njack) '.mat'],'pfinal');

    else
        %save(['ModelCell_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '_random.mat'],'pfinal');
    end
end

disp(['Optimization took ' num2str(toc/60) ' minutes']);

%Plot the results and save the figures


h=pfinal(2:Nlags*Nf+1);


 J=pfinal(Nlags*Nf+2:end); % this is the covariance matrix J
 [V,D] = eig(reshape(J,Nlags*Nf,Nlags*Nf));


figure
subplot(3,3,1)
 imagesc(reshape(h,Nf,Nlags))
 axis xy
 
subplot(3,3,2)
imagesc(reshape(J,Nlags*Nf,Nlags*Nf))
axis xy
axis square

subplot(3,3,3)
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');

subplot(3,3,4)
eig_sorted_1 = V(:,index(1));
imagesc(reshape(eig_sorted_1,Nf,Nlags))
axis xy

subplot(3,3,5)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,Nf,Nlags))
axis xy

subplot(3,3,6)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,Nf,Nlags))
axis xy

subplot(3,3,7)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,Nf,Nlags))
axis xy

subplot(3,3,8)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,Nf,Nlags))
axis xy

subplot(3,3,9)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,Nf,Nlags))
axis xy



figure_name = (['st_' num2str(cellnum) 'your figure name' '_Nlags' num2str(Nlags) '_nfft128_Nf16' '_jack_' num2str(jack) '_of_' num2str(njack)]); % the 6 songs

hgsave(figure_name); 


end
end
