function [pbest, flist, ftestlist] = frprmn(p, func, dfunc, stim, resp, teststim, testresp, order, avgs, Nd, fittype)

ITMAX = 1000;
fp = feval(func, p, stim, resp, order);
xi = feval(dfunc, p, stim, avgs, order);
exitCondition = 0;
g  = -xi;
h  = g;
xi = g;
besttest = 1000;
flist=[];
ftestlist=[];
tally = 0;

% Loop over iterations of minimization
for its=1:ITMAX,
    %disp(['Iteration ' num2str(its)]);
  
    [p, xi, fret] = dlinmin(p, xi, func, dfunc, stim, resp, order, avgs);
    flist(its)=fret;
    if fittype==0
        ftestlist(its)=feval(func, p, teststim, testresp, order);
    end
    
    
%    figure(1)
%    plot(flist)
%    if fittype==0
%        hold on
%        plot(ftestlist,'r')
%        hold off
%    end
%    drawnow

    
    if fittype==0
        
        if ftestlist(its)<besttest*.999999  || its<=2
            besttest = ftestlist(its);
            pbest = p;
            tally=0;
        else
            tally = tally+1;
        end
        
        if tally==10 || its==400
            disp('min of test set found');
            exitCondition = 1;
            break;
        end
        
    else
        J = reshape(p(52:2551),50,50);
        [evecs,evals]=eig(J);
        [EV,inds] = sort((diag(evals)));
        disp(num2str([min(EV) , max(EV)]));
        if its==200
            pbest = p;
            disp('stopping algorithm');
            exitCondition = 1;
            break;
        end
    end
    
    xi = feval(dfunc, p, stim, avgs, order);
    gg = sum(g.^2);
    dgg = sum( (xi + g).*xi );   % This statement for Polak-Ribiere
    % dgg = sum( xi.^2);         % This statement for Fletcher-Reeves
    if gg == 0,            % Unlikely.  If gradient is exactly zero then
        exitCondition = 2;   % we are already done.
        disp('Gradient equal to zero, exiting frprmn.');
        break;
    end
    gam = dgg/gg;
    g = -xi;
    h = g + gam.*h;
    xi = h;
end
if exitCondition == 0,
    disp('Too many iterations in frprmn');
end
