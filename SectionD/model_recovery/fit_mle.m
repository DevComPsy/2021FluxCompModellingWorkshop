function [best_fits] = fit_mle(cfg, fit_function)

%determine likelihood function
if fit_function == 1
    fh = str2func('lik_1LR');
    fit.nparams = 2; %tau, alpha
    fit.lb = [.1, .01];
    fit.ub = [30, 1];
elseif fit_function == 2
    fh = str2func('lik_2LR');
    fit.nparams = 3; %tau, alphapos, alphaneg
    fit.lb = [.1, .01, .01];
    fit.ub = [30, 1, 1];
end

%set number of iterations
niter = 5; 

%get participant data
resp     = cfg.resp; %get responses
ntrl    = length(resp); %get number of trials
rewards = cfg.rewards; % get rewards


%fit model niter times
for iter = 1:niter
    
% choose a random number between the lower and upper bounds to initialize each of the parameters
 fit.init(iter, :) = rand(1,length(fit.lb)).* (fit.ub - fit.lb) + fit.lb;
 
%fit model
[res, nll] =  ...
    fmincon(@(x) fh(resp, rewards, x),...
    squeeze(fit.init(iter,:)),[],[],[],[],fit.lb, fit.ub,[],...
    optimset('maxfunevals',10000,'maxiter',2000, 'Display', 'off'));

% Store results for each subject
    fit.result.lik(iter) = nll;
    fit.result.AIC(iter) = 2*nll + fit.nparams*2;
    fit.result.BIC(iter) = 2*nll + fit.nparams*log(ntrl);
    fit.result.params(iter, :) = res;

 %Find iteration with minimum neg log lik
    [a, b] = min(fit.result.lik);
        
    %save best fitting parameters
    best_fits(:) = [fit.result.lik(b), fit.result.AIC(b), fit.result.BIC(b), fit.result.params(b, :)];
end
end

          