function pfinal = MNEfit(stim, resp, teststim, testresp, celltype, cellnum, jack, order, savepath)


% Calculate constraints
[Nsamples,Ndim] = size(stim);
avg = (stim'*resp)/Nsamples;  % Ndim x Nrep
avg = mean(avg,2);  % Ndim x 1
avgsqrd = stim'*(repmat(resp,[1,Ndim]).*stim)/Nsamples;  % Ndim x Ndim
avgsqrd = reshape(avgsqrd,[Ndim^2,1]);
psp = mean(mean(resp));   % spike probability
avgs = [psp;avg];
if order>1
    avgs = [avgs;avgsqrd];
end



% Create a starting point
pstart = log(1/avgs(1)-1);
%pstart = 0;
pstart(2:Ndim+1) = .000001*(2*rand([1,Ndim])-1);
if order>1
    temp = .0000005*(2*rand([Ndim,Ndim])-1);
    pstart(Ndim+2:length(pstart)+Ndim^2) = reshape((temp+temp'),[1,Ndim^2]);
    clear temp;
end
ftol=1;



% Run conjugate gradient algorithm
[pfinal, flist, ftestlist] = frprmnRGC(pstart, ftol, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, savepath, Nd);



save([savepath celltype '_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '.mat'],'pfinal', '-v6');
%save([savepath 'flist_' celltype '_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '.mat'],'flist', '-v6');
%save([savepath 'ftestlist_' celltype '_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '.mat'],'ftestlist', '-v6');


