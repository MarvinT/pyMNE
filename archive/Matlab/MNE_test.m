function MNE_test()
    load('example_MNE.mat') 
    %{
    I only include a minimal amount of data so that it runs fast
    contains:
        example['stimuli'] = stimuli
        example['spikes'] = spikes
        example['pstart'] = pstart
        example['v'] = v
        example['ll_pstart'] = logLoss(pstart)
        example['dll_pstart'] = dlogLoss(pstart)
        example['pstart2'] = pstart2
        example['v'] = v
        example['ll_pstart2'] = logLoss2(pstart2)
        example['dll_pstart2'] = dlogLoss2(pstart2)
    %}
    order = 1;
    Nd = 0;
    fittype = 0;

    [Nsamples,Ndim] = size(stimuli);
    psp = mean(mean(spikes));   % spike probability
    avg = (stimuli'*spikes)/Nsamples;  % Ndim x Nrep
    avg = mean(avg,2);  % Ndim x 1 I dont think this does anything
    avgs = [psp;avg];
    if order>1
        avgsqrd = stimuli'*(repmat(spikes,[1,Ndim]).*stimuli)/Nsamples;  % Ndim x Ndim
        avgsqrd = reshape(avgsqrd,[Ndim^2,1]);
        avgs = [avgs;avgsqrd];
    end
    
    fstart = logloss(pstart, stimuli, spikes, order);
    dfstart = dlogloss(pstart, stimuli, avgs, order);
    
    display('difference in evaluated log loss')
    display( fstart - ll_pstart )
    display('difference in evaluated d log loss')
    display( max(dfstart - dll_pstart) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matlab code from fitzgerald
% this should all be correct, I'm just including it for completeness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = logloss(p, stim, resp, order)
    [Nsamples,Ndim] = size(stim);

    ptemp = p(2:Ndim+1);
    if order>1
        J = reshape(p(Ndim+2:Ndim+1+Ndim^2),[Ndim,Ndim])';
    end

    if order==1
        f1 = 1+exp(p(1)+stim*ptemp');
        f0 = 1+exp(-p(1)-stim*ptemp');
    else
        f1 = 1+exp(p(1)+stim*ptemp'+sum(stim.*(stim*J),2));
        f0 = 1+exp(-p(1)-stim*ptemp'-sum(stim.*(stim*J),2));
    end

    NspikesperStim = resp;  % Nsamples x 1
    F1 = NspikesperStim.*log(f1);
    F0 = (1-NspikesperStim).*log(f0);
    F1(isnan(F1)) = 0;
    F0(isnan(F0)) = 0;
    f = mean(F0+F1);
end

function df = dlogloss(p, stim, avgs, order)
    [Nsamples,Ndim] = size(stim);

    ptemp = p(2:Ndim+1);
    if order>1
        J = reshape(p(Ndim+2:Ndim+1+Ndim^2),[Ndim,Ndim]);
    end

    if order==1
        pSpike = 1./(1+exp(p(1)+stim*ptemp'));  % Nsamples x 1
        averages = mean(pSpike);
        averages(2:Ndim+1,1) = stim'*pSpike/Nsamples;
    elseif order==2
        pSpike = 1./(1+exp(p(1)+stim*ptemp'+sum(stim.*(stim*J),2)));  % Nsamples x 1
        averages = mean(pSpike);
        averages(2:Ndim+1,1) = stim'*pSpike./Nsamples;
        temp = stim'*(repmat(pSpike,[1,Ndim]).*stim)./Nsamples;  % Ndim x Ndim
        temp = reshape(temp,[Ndim^2,1]);
        averages(Ndim+2:Ndim+1+Ndim^2) = temp;    
    end

    df = (avgs - averages)';
end