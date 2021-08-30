%% KN
% Likelihood function for one LR model

function [lik] = null_lik(blocks, choices, rewards, x)

% Initialize log likelihood at 0
lik = 0;

% Loop through trials
for trial = 1:length(choices)
    
    choice_probs = [.5 .5];
    
    % Determine the choice the participant actually made on this trial
    trial_choice = choices(trial);
    
    %Determine the probability that the participant made the choice they
    %made
    lik_choice = choice_probs(trial_choice);
    
    %update log likelihood
    lik = lik + log(lik_choice);
    
 
end

% Put priors on parameters 
%lik= lik+log(pdf('beta', alpha, 1.1,1.1));
%lik= lik+log(pdf('gam', beta, 4.82, 0.88));

%flip sign of loglikelihood (which is negative, and we want it to be as close to 0 as possible) so we can enter it into fmincon, which searches for minimum, rather than maximum values
lik = -lik;

