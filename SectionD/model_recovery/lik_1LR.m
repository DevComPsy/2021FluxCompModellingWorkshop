%% KN
% Likelihood function for one LR model

function [lik] = lik_1lr(resp, rewards, x)

ntrl = length(resp);
tau = x(1);
alpha = x(2);

% Initialize q values
q = nan(ntrl, 2);
q(1,:) = 0.5; %set at .5 for first trial
vm =  [.5, .5]; %set cached value for each bandit;

% Loop through trials
for itrl = 2:ntrl
    
    %get reward and update cached values
    r = resp(itrl-1);
    vm(r) = rewards(r, itrl-1);
    vm(3-r) = 1-vm(r);
    
    %update q values
    q(itrl,r)   = q(itrl-1,r)+alpha*(vm(r)-q(itrl-1,r));
    q(itrl,3-r) = q(itrl-1,3-r);
    
end

%get log of choice probabilities on every trial 
probs = 1./(1+exp(-(q(:,1)-q(:,2))/tau)); 
log_choice_probs = log([probs(resp == 1); 1-probs(resp == 2)]);

%compute log likelihood
lik = sum(log_choice_probs);

%flip sign of loglikelihood (which is negative, and we want it to be as close to 0 as possible) so we can enter it into fmincon, which searches for minimum, rather than maximum values
lik = -lik;

