% Script to simulate and test recoverability of three models:
% 1 LR
% 2 LR
% Null model


% Kate Nussenbaum (adapted from Vasilisa Skvortsova)
% 08/2021 
% FLUX Workshop on Computational modeling

clear all; 
close all; 
clc; 

%% Task presets
% Preset the task: 2-armed bandit task with binary rewards and probability
% of reward 0.8
% Simulate N = 1000 subjects
cfg_sim.nsims      = 1000; %number of parameter combinations
cfg_sim.nsub      = 1; %number of participants with each parameter combination
cfg_sim.p         = 0.8; % 80 20 
cfg_sim.nblck     = 1; % assuming 1 block 
cfg_sim.ntrl      = 150; % assume we have 2 blocks of 60 trials each with 2 reversals within each block (for simplicity in the middle of the block).

heatmap_name = 'model_recoverability_150.png';

%% PART 1: Simulate  data from all three models 
% 1 LR data
for idx = 1:cfg_sim.nsims
    fprintf('Running simulation %d\n',idx)
    cfg_sim.tau = unifrnd(0.1, 1.5);
    cfg_sim.alphapos = unifrnd(0.05, 1); %sample positive learning rate from between .05 and 1
    cfg_sim.alphaneg = unifrnd(0.05, 1);  % neg isn't used in 1 LR model
    sim_out  = simulate_1LR_data(cfg_sim);
    sim_data_1LR.(sprintf('sim%d',idx)) = sim_out; 
end

% 2 LR data
for idx = 1:cfg_sim.nsims
    fprintf('Running simulation %d\n',idx)
    cfg_sim.tau = unifrnd(0.1, 1.5);
    cfg_sim.alphapos = unifrnd(0.4, 0.5); %sample positive learning rate from between .4 and .5
    cfg_sim.alphaneg = unifrnd(0.05, .1);  %sample negative learning rate from between .05 and .1
    sim_out  = simulate_2LR_data(cfg_sim);
    sim_data_2LR.(sprintf('sim%d',idx)) = sim_out; 
end

% null data
for idx = 1:cfg_sim.nsims
    fprintf('Running simulation %d\n',idx)
    cfg_sim.tau = unifrnd(0.1, 1.5);
    cfg_sim.alphapos = unifrnd(0.05, 1); %sample positive learning rate from between .05 and 1
    cfg_sim.alphaneg = unifrnd(0.05, 1);  %sample negative learning rate from between .05 and 1
    sim_out  = simulate_null_data(cfg_sim);
    sim_data_null.(sprintf('sim%d',idx)) = sim_out; 
end


%% PART 2: Fit each model to each dataset
clc; 

cfg_fit.nsims = cfg_sim.nsims;
fit_data = []; 
out_1LR = [];

for idx = 1:cfg_fit.nsims

        %fit models to 1 LR data
        fprintf('Fitting to one LR data %d\n',idx)
        cfg_fit.resp    = sim_data_1LR.(sprintf('sim%d',idx)).resp; %get responses
        cfg_fit.rewards = sim_data_1LR.(sprintf('sim%d',idx)).rewards; %get rewards
        
        %fit models to data and store results
        out_1LRdata_1LR(idx, :) = fit_mle(cfg_fit, 1); 
        out_1LRdata_2LR(idx, :) = fit_mle(cfg_fit, 2);
        
        %fit models to 2 LR data
        fprintf('Fitting to two LR data %d\n',idx)
        cfg_fit.resp    = sim_data_2LR.(sprintf('sim%d',idx)).resp; %get responses
        cfg_fit.rewards = sim_data_2LR.(sprintf('sim%d',idx)).rewards; %get rewards
        
        %fit models to data and store results
        out_2LRdata_1LR(idx, :) = fit_mle(cfg_fit, 1); 
        out_2LRdata_2LR(idx, :) = fit_mle(cfg_fit, 2); 
        
        %fit models to null data
        fprintf('Fitting to null data %d\n',idx)
        cfg_fit.resp    = sim_data_null.(sprintf('sim%d',idx)).resp; %get responses
        cfg_fit.rewards = sim_data_null.(sprintf('sim%d',idx)).rewards; %get rewards
        
        %fit models to data and store results
        out_nulldata_1LR(idx, :) = fit_mle(cfg_fit, 1); 
        out_nulldata_2LR(idx, :) = fit_mle(cfg_fit, 2); 
end

%% PART 3: Plot recoverability
make_plots;


%% === Additional functions ====== %%
function out = simulate_1LR_data(cfg)

if ~isfield(cfg,'p')
    fprintf('No reward probability specified, assume probability of reward is 0.75\n')
    p = 0.75;
else
    p = cfg.p;
end


if ~isfield(cfg,'tau')
    fprintf('No value for softmax, assuming tau = 0.1\n')
    tau = 0.1;
else
    tau  = cfg.tau;
end


if ~isfield(cfg,'alphapos')
    fprintf('No value for postiive learning rate, assuming learnng rate chosen = 0.5\n')
    alphapos = 0.5;
else
    alphapos  = cfg.alphapos;
end


if ~isfield(cfg,'alphaneg')
    fprintf('No value for negative learning rate, assuming same as positive learning rate\n')
    alphaneg = cfg.alphapos;
else
    alphaneg  = cfg.alphaneg;
end

if ~isfield(cfg,'nsub')
    fprintf('Assuming 1 subject\n')
    nsub  = 1;
else
    nsub  = cfg.nsub;
end


if ~isfield(cfg,'ntrl')
    fprintf('Assuming 120 trials\n')
    
    ntrl  = 120;
else
    ntrl  = cfg.ntrl;
end


if ~isfield(cfg,'nblck')
    fprintf('Assuming 2 blocks\n')
    
    nblck  = 2;
else
    nblck  = cfg.nblck;
end


ntrl_blck = ntrl/nblck;



for i_s = 1:nsub
    
    fprintf('Simulating subject %d\n',i_s)
    % create binary rewards with reward probability p
    rew1 = [];
    rew2 = [];
    
    for i_b = 1:nblck
        
        r1                  = zeros(1,ntrl);
        r2                  = zeros(1,ntrl); 
        
        r1(1:floor(ntrl*p))     = 1; 
        r1(1:floor(ntrl*(1-p))) = 1; 
        
        r1 = r1(randperm(length(r1))); 
        r2 = r2(randperm(length(r2))); 
        
        if mod(i_b,2) == 1
            rew1 = cat(2,rew1,r1);
            rew2 = cat(2,rew2,r2);
            
        else
            rew1 = cat(2,rew1,r2);
            rew2 = cat(2,rew2,r1);
        end
    end
    
    rewards = cat(1,rew1,rew2);
    
    resp    = get_simulate_1LR(rewards,alphapos,alphaneg,tau); 
    
    out(i_s).resp    = resp; 
    out(i_s).rewards = rewards; 
    
    
%     figure; 
%     plot(out.rewards(1,:),'or'); hold on; plot(out.rewards(2,:),'ob')
%     plot(out.resp-1,'ok');
     
    
end
end



function out = simulate_2LR_data(cfg)

if ~isfield(cfg,'p')
    fprintf('No reward probability specified, assume probability of reward is 0.75\n')
    p = 0.75;
else
    p = cfg.p;
end


if ~isfield(cfg,'tau')
    fprintf('No value for softmax, assuming tau = 0.1\n')
    tau = 0.1;
else
    tau  = cfg.tau;
end


if ~isfield(cfg,'alphapos')
    fprintf('No value for postiive learning rate, assuming learnng rate chosen = 0.5\n')
    alphapos = 0.5;
else
    alphapos  = cfg.alphapos;
end


if ~isfield(cfg,'alphaneg')
    fprintf('No value for negative learning rate, assuming same as positive learning rate\n')
    alphaneg = cfg.alphapos;
else
    alphaneg  = cfg.alphaneg;
end

if ~isfield(cfg,'nsub')
    fprintf('Assuming 1 subject\n')
    nsub  = 1;
else
    nsub  = cfg.nsub;
end


if ~isfield(cfg,'ntrl')
    fprintf('Assuming 120 trials\n')
    
    ntrl  = 120;
else
    ntrl  = cfg.ntrl;
end


if ~isfield(cfg,'nblck')
    fprintf('Assuming 2 blocks\n')
    
    nblck  = 2;
else
    nblck  = cfg.nblck;
end


ntrl_blck = ntrl/nblck;



for i_s = 1:nsub
    
    fprintf('Simulating subject %d\n',i_s)
    % create binary rewards with reward probability p
    rew1 = [];
    rew2 = [];
    
    for i_b = 1:nblck
        
        r1                  = zeros(1,ntrl);
        r2                  = zeros(1,ntrl); 
        
        r1(1:floor(ntrl*p))     = 1; 
        r1(1:floor(ntrl*(1-p))) = 1; 
        
        r1 = r1(randperm(length(r1))); 
        r2 = r2(randperm(length(r2))); 
        
        if mod(i_b,2) == 1
            rew1 = cat(2,rew1,r1);
            rew2 = cat(2,rew2,r2);
            
        else
            rew1 = cat(2,rew1,r2);
            rew2 = cat(2,rew2,r1);
        end
    end
    
    rewards = cat(1,rew1,rew2);
    
    resp    = get_simulate_2LR(rewards,alphapos,alphaneg,tau); 
    
    out(i_s).resp    = resp; 
    out(i_s).rewards = rewards; 
    
    
%     figure; 
%     plot(out.rewards(1,:),'or'); hold on; plot(out.rewards(2,:),'ob')
%     plot(out.resp-1,'ok');
     
    
end
end



function out = simulate_null_data(cfg)

if ~isfield(cfg,'p')
    fprintf('No reward probability specified, assume probability of reward is 0.75\n')
    p = 0.75;
else
    p = cfg.p;
end


if ~isfield(cfg,'tau')
    fprintf('No value for softmax, assuming tau = 0.1\n')
    tau = 0.1;
else
    tau  = cfg.tau;
end


if ~isfield(cfg,'alphapos')
    fprintf('No value for postiive learning rate, assuming learnng rate chosen = 0.5\n')
    alphapos = 0.5;
else
    alphapos  = cfg.alphapos;
end


if ~isfield(cfg,'alphaneg')
    fprintf('No value for negative learning rate, assuming same as positive learning rate\n')
    alphaneg = cfg.alphapos;
else
    alphaneg  = cfg.alphaneg;
end

if ~isfield(cfg,'nsub')
    fprintf('Assuming 1 subject\n')
    nsub  = 1;
else
    nsub  = cfg.nsub;
end


if ~isfield(cfg,'ntrl')
    fprintf('Assuming 120 trials\n')
    
    ntrl  = 120;
else
    ntrl  = cfg.ntrl;
end


if ~isfield(cfg,'nblck')
    fprintf('Assuming 2 blocks\n')
    
    nblck  = 2;
else
    nblck  = cfg.nblck;
end


ntrl_blck = ntrl/nblck;



for i_s = 1:nsub
    
    fprintf('Simulating subject %d\n',i_s)
    % create binary rewards with reward probability p
    rew1 = [];
    rew2 = [];
    
    for i_b = 1:nblck
        
        r1                  = zeros(1,ntrl);
        r2                  = zeros(1,ntrl); 
        
        r1(1:floor(ntrl*p))     = 1; 
        r1(1:floor(ntrl*(1-p))) = 1; 
        
        r1 = r1(randperm(length(r1))); 
        r2 = r2(randperm(length(r2))); 
        
        if mod(i_b,2) == 1
            rew1 = cat(2,rew1,r1);
            rew2 = cat(2,rew2,r2);
            
        else
            rew1 = cat(2,rew1,r2);
            rew2 = cat(2,rew2,r1);
        end
    end
    
    rewards = cat(1,rew1,rew2);
    
    resp    = get_simulate_null(rewards,alphapos,alphaneg,tau); 
    
    out(i_s).resp    = resp; 
    out(i_s).rewards = rewards; 
    
    
%     figure; 
%     plot(out.rewards(1,:),'or'); hold on; plot(out.rewards(2,:),'ob')
%     plot(out.resp-1,'ok');
     
    
end
end



% simulate agent with specific learning parameters

% 1 LR model
function [resp] = get_simulate_1LR(rewards,alphapos,alphaneg,tau)

    % compute Q-values
    alpha = alphapos;
    ntrl    = size(rewards,2); 
    q       = nan(ntrl,2);
    resp    = nan(ntrl,1); 
    
    q(1,:) = 0.5;
    vm     = [0.5,0.5]; % cached value for each bandit
    
    % first choice is random 
    resp(1) = double(rand>0.5)+1; % between 1 and 2 
    
    for itrl = 2:ntrl
        
        r = resp(itrl-1); %get response from previous trial;
        vm(r)   = rewards(r,itrl-1); % get reward associated with response
        
        % update the Q-values
        q(itrl,r)   = q(itrl-1,r)+alpha*(vm(r)-q(itrl-1,r)); %update chosen option
        q(itrl,3-r) = q(itrl-1,3-r); % don't update unchosen option
        
        
        deltaQ = q(itrl,1) - q(itrl,2);    % difference in Q values between options 1 and option 2 
        proba  = 1./(1+exp(-deltaQ/tau));  % softmax for 2 options 
        
        resp(itrl) = (1-sampleFromArbitraryP([proba,1-proba],[1,0]',1))+1; %get response
        
    end
   
end


% 2 LR model
function [resp] = get_simulate_2LR(rewards,alphapos,alphaneg,tau)

    % compute Q-values
    ntrl    = size(rewards,2); 
    q       = nan(ntrl,2);
    resp    = nan(ntrl,1); 
    
    q(1,:) = 0.5;
    vm     = [0.5,0.5]; % cached value for each bandit
    
    % first choice is random 
    resp(1) = double(rand>0.5)+1; % between 1 and 2 
    
    for itrl = 2:ntrl
        
        r = resp(itrl-1);
        vm(r)   = rewards(r,itrl-1); % update cached value for chosen bandit
        %vm(3-r) = 1 - vm(r); % assume anticorrelation
        
        % update the Q-values for chosen option
        if (vm(r)-q(itrl-1,r))> 0
            q(itrl,r)   = q(itrl-1,r)+alphapos*(vm(r)-q(itrl-1,r));
        else
            q(itrl,r)   = q(itrl-1,r)+alphaneg*(vm(r)-q(itrl-1,r));
        end
        
        % don't update unchosen option
        q(itrl,3-r) = q(itrl-1,3-r); 
        
        deltaQ = q(itrl,1) - q(itrl,2);    % difference in Q values between options 1 and option 2 
        proba  = 1./(1+exp(-deltaQ/tau));  % softmax for 2 options 
        
        resp(itrl) = (1-sampleFromArbitraryP([proba,1-proba],[1,0]',1))+1; 
        
    end
   
end


% Null model
function [resp] = get_simulate_null(rewards,alphapos,alphaneg,tau)

    ntrl    = size(rewards,2); 
    resp    = nan(ntrl,1); 
    for itrl = 1:ntrl
        
        % all choices are random
        resp(itrl) = double(rand>0.5)+1; % between 1 and 2 
      
    end
   
end











