%% Generate simulated choice data %%
% Kate Nussenbaum - katenuss@nyu.edu
% Flux Computational Modeling Workshop
% Last edited: 8/30/21

%clear everything
clear;

%reset random stream
rng('shuffle');

%load path to models
addpath('models');

%% VARIABLES TO MODIFY %%
% Data file name
save_filename = ['simulated_data/data_', date];

% Number of subjects to simulate
num_subs = 1000;

% Task structure
task_struct.reward_probs = [.8 .2]; %reward probabilities
task_struct.num_blocks = 3; %number of blocks (assumes value estimates are re-set at beginning of each block)
task_struct.num_block_trials = 50; %number of trials per block

%Models to simulate
models = {'1LR', '2LR', 'decay', 'null'};

% Note: Can also modify parameter distributions in code below

%% GENERATE SIMULATED DATA %%
%----------------------------------%
% Loop through models and subjects %
%----------------------------------%
for m = 1:length(models)
    model_to_simulate = models{m};
    
    %print message about which subject is being fit
    fprintf('Simulating model %d out of %d...\n', m, length(models))
    
    % draw parameters for each model
    if strcmp(model_to_simulate, '1LR')
        tau = rand(1, num_subs); %between 0 and 1
        alpha = rand(1, num_subs); %between 0 and 1
        params = [tau', alpha'];
        param_names = {'beta', 'alpha'};
        function_name = 'sim_oneLR';
    elseif strcmp(model_to_simulate, '2LR') %for simulation purposes, make the learning rates different from each other
        tau = rand(1, num_subs);  %between 0 and 1
        alpha_pos = .2 * rand(1, num_subs) + .3;  %between .3 and 5
        alpha_neg = .05 * rand(1, num_subs) + .05; %between 0.05 and .1
        params = [tau', alpha_pos', alpha_neg'];
        param_names = {'beta', 'alpha_pos', 'alpha_neg'};
        function_name = 'sim_twoLR';
    elseif strcmp(model_to_simulate, 'null')
        params = [];
        param_names = {};
        function_name = 'sim_null';
    elseif strcmp(model_to_simulate, 'decay')
        tau = rand(1, num_subs);  %between 0 and 1
        alpha_init = rand(1, num_subs);  %between 0 and 1
        eta = rand(1, num_subs); %between 0 and 1
        params = [tau', alpha_init', eta'];
        param_names = {'beta', 'alpha_init', 'eta'};
        function_name = 'sim_decay';
    end
    
    % convert function name into function
    fh = str2func(function_name);
    
    %loop through subjects and generate data
    for s = 1:num_subs
        if isempty(params)
            [blocks, choices, rewards] = fh(task_struct);
            model_data(s).params = [];
        else
            [blocks, choices, rewards] = fh(task_struct, params(s, :));
            model_data(s).params = params(s, :);
            
        end
        
        %for each subject, save the generating parameters, choices, and
        %rewards
        model_data(s).choices = choices;
        model_data(s).rewards = rewards;
        model_data(s).blocks = blocks;
    end
    
    %for each model, save the simulated subject data, name of the model, 
    %number of parameters, and parameter names
    sim_data(m).sub_data = model_data;
    sim_data(m).function = function_name;
    sim_data(m).n_params = size(params, 2);
    sim_data(m).param_names = param_names;
    
end

%---------------------%
% Save simulated data %
%---------------------%
save(save_filename, 'sim_data');
