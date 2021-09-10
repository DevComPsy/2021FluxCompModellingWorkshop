%% Fit models to simulated choice data %%
% Kate Nussenbaum - katenuss@nyu.edu
% Flux Computational Modeling Workshop
% Last edited: 8/30/21

% clear everything
clear;

%load path to likelihood functions
addpath('lik_funs');

%% VARIABLES TO MODIFY %%
%Data file name
save_filename = ['simulated_data/200_trials_fits_', date];

% Data to load
load('simulated_data/200_trials_data_10-Sep-2021');
data = sim_data;

% How many iterations to run per participant
niter = 5;

% Models to fit
models = {'1LR', 'decay', 'null'};

%determine the number of subjects
n_subjects = length(data(1).sub_data);

%% FIT MODELS TO DATA %%

%---------------------------------------------%
% Loop through models, datasets, and subjects %
%---------------------------------------------%
for d = 1:size(sim_data, 2)
    
    %get data for this dataset
    dataset_data = data(d).sub_data;
    
    %print message about which dataset is being fit
    fprintf('Fitting dataset %d out of %d...\n', d, size(sim_data, 2))
    
    for m = 1:length(models)
        model_to_fit = models{m};
        
        %print message about which subject is being fit
        fprintf('Fitting model %d out of %d...\n', m, length(models))
        
        % model specific info
        if strcmp(model_to_fit, '1LR')
            n_params = 2; %tau, alpha
            lb = [.1, 1e-6]; %lower bounds of the parameters
            ub = [3, 1]; %upper bounds of the parameters
            function_name = 'one_LR_lik';
        elseif strcmp(model_to_fit, '2LR')
            n_params = 3; %tau, alpha_pos, alpha_neg
            lb = [.1, 1e-6, 1e-6];
            ub = [3, 1, 1];
            function_name = 'two_LR_lik';
        elseif strcmp(model_to_fit, 'decay')
            n_params = 3; %tau, alpha_init, eta
            lb = [.1, 1e-6, 1e-6];
            ub = [3, 1, 1];
            function_name = 'decay_lik';
        elseif strcmp(model_to_fit, 'null')
            n_params = 0; %no parameters
            function_name = 'null_lik';
        end
        
        % convert function name to function
        fh = str2func(function_name);
        
        % generate matrices to save data
        [logpost, negloglik, AIC, BIC] = deal(nan(n_subjects, 1));
        [params] = nan(n_subjects, n_params);
        
        
        %loop through subjects
        parfor s = 1:n_subjects
            sub_data = dataset_data(s);
            
            fprintf('Fitting subject %d out of %d...\n', s, n_subjects) %print message saying which subject is being fit
            
            %if null model, compute negative log likelihood directly
            if strcmp(model_to_fit, 'null')
                logpost(s) = NaN;
                params(s, :) = NaN;
                negloglik(s) = -1*sum(log(.5*ones(length(sub_data.choices), 1)));
                AIC(s) = 2*negloglik(s);
                BIC(s) = 2*negloglik(s);
            else %otherwise fit models
                for iter = 1:niter  % run niter times from random initial conditions, to get best fit
                    % choosing a random number between the lower and upper bounds
                    % (defined above) to initialize each of the parameters
                    starting_points = rand(1,length(lb)).* (ub - lb) + lb; % random initialization
                    
                    % Run fmincon (unless null model)
                    [res,nlp] = ...
                        fmincon(@(x) fh(sub_data.blocks, sub_data.choices, sub_data.rewards, x),...
                        starting_points,[],[],[],[],lb, ub,[],...
                        optimset('maxfunevals',10000,'maxiter',2000, 'Display', 'off'));
                    
                    %flip sign to get log posterior (if priors are in models, if no priors, this will be the log likelihood)
                    logp = -1 * nlp;
                    
                    %store results if minimum is found
                    if iter == 1 || logpost(s) < logp
                        logpost(s) = logp;
                        params(s, :) = res;
                        negloglik(s) = fh(sub_data.blocks, sub_data.choices, sub_data.rewards, res); %fit model w/ 'wining' parameters to get the negative log likelihood
                        AIC(s) = 2*negloglik(s) + 2*length(res);
                        BIC(s) = 2*negloglik(s) + length(res)*log(length(sub_data.choices));
                    end
                    
                end
            end
        end
        
        results.logpost = logpost;
        results.params = params;
        results.negloglik = negloglik;
        results.AIC = AIC;
        results.BIC = BIC;
        
        %save for each model
        model_fits(d, m).results = results;
        model_fits(d, m).fit_model = model_to_fit;
        model_fits(d, m).sim_model = data(d).function;
        
    end
end


%----------------------------%
% Save model-fitting results %
%----------------------------%
save(save_filename, 'model_fits');



