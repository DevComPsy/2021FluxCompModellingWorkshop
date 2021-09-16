%-------------------------------------------------------------------------%
%    SIMULATION OF THE SIMPLE Q-LEARNING MODEL FOR INSTRUMENTAL CONDITIONING
%
%-------------------------------------------------------------------------%
clear all
close all


rng('shuffle')
%% -------------------------- INPUTS  -------------------------------------

% PARAMETERS
alpha      = 0.2;
inv_temp   = 3;

% NUMBER OF TRIALS
ntrials = 24;
nsess  = 50;


%% ------------------------- PARAMETERS -----------------------------------

% INITIAL VALUE
Q0  = [0.5 0.5]; % initial value


%% ------------------------- Q SIMULATION --------------------------------

% initialise value and prediction error vectors to store the data
Qt  = nan(ntrials+1, 2,nsess);
PE  = nan(ntrials, nsess);
ch  = nan(ntrials, nsess);
PA  = nan(ntrials, nsess);
r   = nan(ntrials, nsess);


for ksess = 1:nsess
   
    % gemerate outcomes & initialise Q vector
   O = [rand(ntrials,1)<0.2, rand(ntrials,1)<0.8];
   Qt(1,:,ksess)  = Q0;
    
    % value simulation
    for t = 1:ntrials
        
        PA(t,ksess)   = 1./(1+exp(-inv_temp.*(Qt(t,2,ksess)-Qt(t,1,ksess))));
        ch(t,ksess)   = 1 + double(rand()<PA(t,ksess));
        
        % outcome
        r(t,ksess) = O(t,ch(t,ksess));
        
        % compute prediction error
        PE(t,ksess) = r(t,ksess) - Qt(t,ch(t,ksess),ksess);
        
        % update value
        Qt(t+1,ch(t,ksess),ksess) = Qt(t,ch(t,ksess),ksess) + alpha.*PE(t,ksess);
        Qt(t+1,3-ch(t,ksess),ksess) = Qt(t,3-ch(t,ksess),ksess);
        
    end
end
%% ------------------------- PLOTS  ---------------------------------------
figure

% display the simulation results for predicted value
subplot(2, 1, 1)
hold on
plot(squeeze(mean(Qt(:,1,:),3)), '-b');
plot(squeeze(mean(Qt(:,2,:),3)), '-r');
xlim([0, ntrials]);
ylim([0, 1.1]);
xlabel('trials');
ylabel('Q-vvalues');
legend('A','B');

% display the simulation results for prediction error
subplot(2, 1, 2);
hold on
plot(squeeze(mean(PA,2)), '-k');
plot(squeeze(mean(ch-1,2)), 'o');
xlim([0, ntrials]);
ylim([-0.1, 1.1]);
xlabel('trials');
ylabel('Proba(choose A)');
legend('proba','choice');

alphaSim    = alpha;
betaSim     = inv_temp;
save('sim_data_sess','ch','r','nsess','alphaSim','betaSim')
