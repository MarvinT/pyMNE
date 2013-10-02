function p = RGC_LGN_RunMNEfit(celltype, cellnum, jack, order, stimtype, fittype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function loads the stimulus and spikes, breaks them into the %
% appropriate jackknife and starts the optimization                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Running ' celltype num2str(cellnum) ' MNE' num2str(order) ' jackknife ' num2str(jack)]);

spikedir = ['./Data/Spikes/' celltype '/'];
stimdir  = ['./Data/Stimuli/RGC_LGN/'];
savepath = ['./Results/RGC_LGN/'];

%%%%%%%%%%%%%%%%%
% load stimulus %
%%%%%%%%%%%%%%%%%
id=fopen([stimdir 'Uniques/S' num2str(cellnum) '.raw'],'rb');
[stim1]=fread(id,'double');
M1=length(stim1);
fclose(id);
id=fopen([stimdir 'Repeats/S' num2str(cellnum) '.raw'],'rb');
[stim2]=fread(id,'double');
M2=length(stim2);
fclose(id);
stimunits = 1;
% either z-scored or normalized to 0-1 range
if stimunits==1
    stim1 = (stim1-mean([stim1;stim2]))/std([stim1;stim2]);
    stim2 = (stim2-mean([stim1;stim2]))/std([stim1;stim2]);
else
    stim1 = (stim1)/max([stim1;stim2]);
    stim2 = (stim2)/max([stim1;stim2]);
end
% reshape stimulus into 50-D column for each time point
stim = zeros(M1-49,50);
for j=1:M1-49
    stim(j,:)=stim1(j:j+49)';
end
[M,Ndim] = size(stim);

%%%%%%%%%%%%%%%%%%
% load responses %
%%%%%%%%%%%%%%%%%%
fid = fopen([spikedir 'Uniques/r' num2str(cellnum) '.isk']);
resp = fread(fid);
fclose(fid);
resp(resp>1)=1;     % enforce binary responses
resp = resp(50:M1); % cut early responses b/c no stimulus for those

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to 2-D MID space if desired %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimtype==1
    MIDdir = ['./Results/RGC_LGN_MID/' celltype '/Avg_Filters/'];
    load([MIDdir celltype '_' num2str(cellnum) '_MID1'])
    load([MIDdir celltype '_' num2str(cellnum) '_MID2'])
    u1 = stim*mid1;
    u2 = stim*mid2;
    clear stim;
    stim = [u1,u2];
    clear u1 u2;
    [M,Ndim] = size(stim);
    Mstim=max(max(stim));
    stim = stim/Mstim;
    Nd=sqrt(2);
    
    % convert the repeated stimulus into MID space and save for future use
    saverepstim = 0;
    if saverepstim==1
        stimrep = zeros(M2-49,50);
        for j=1:M2-49
            stimrep(j,:)=stim2(j:j+49)';
        end
        u1 = stimrep*mid1;
        u2 = stimrep*mid2;
        clear stimrep;
        stimrep = [u1,u2]./Mstim;
        clear u1 u2;
        save([stimdir 'Repeats/s2_' num2str(cellnum) '.mat'],'stimrep');
    end
else
    Nd=sqrt(50);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% break into training and test sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fittype==0
    Ntest = floor(M/4);
    teststim = stim((jack-1)*Ntest+1:jack*Ntest,:);
    stim((jack-1)*Ntest+1:jack*Ntest,:) = [];
    testresp = resp((jack-1)*Ntest+1:jack*Ntest);
    resp((jack-1)*Ntest+1:jack*Ntest) = [];
end
% shuffle the responses if fitting randomly
if fittype==1
    Ntest = floor(M/8);
    teststim = stim((jack-1)*Ntest+1:jack*Ntest,:);
    stim((jack-1)*Ntest+1:jack*Ntest,:) = [];
    testresp = resp((jack-1)*Ntest+1:jack*Ntest);
    clear resp stim
    stim = teststim;
    resp = testresp;
    temp = randperm(length(resp));
    resp = resp(temp);
    temp = randperm(length(testresp));
    testresp = testresp(temp);
end


% start the fitting
tic
p = MNEfit(stim, resp, teststim, testresp, celltype, cellnum, jack, order, savepath, Nd, fittype);
disp(['Optimization took ' num2str(toc/60) ' minutes']);