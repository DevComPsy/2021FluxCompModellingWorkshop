%-------------------------------------------------------------------------%
%    SIMULATION OF THE SIMPLE Q-LEARNING MODEL FOR INSTRUMENTAL CONDITIONING
%
%-------------------------------------------------------------------------%


%% -------------------------- INPUTS  -------------------------------------

% PARAMETERS
alpha      = 0.3;
inv_temp   = 1.5;

% NUMBER OF TRIALS
ntrials = 24;


%% ------------------------- PARAMETERS -----------------------------------

% INITIAL VALUE
Q0  = [0.5 0.5]; % initial value

% REWARD HISTORY
RA = rand(ntrials,1)<0.8;
RB = rand(ntrials,1)<0.2;
O = [RB, RA];


%% ------------------------- Q SIMULATION --------------------------------

% initialise value and prediction error vectors to store the data
Qt  = nan(ntrials+1, 2);    % store the Q-values for option B (col 1) and option A (col 2)
PE  = nan(ntrials, 1);      % store the Prediction error
ch  = nan(ntrials, 1);      % store the choice : 1 = B; 2 = A
PA  = nan(ntrials, 1);      % store the modelled proba. of choosing A
r     = nan(ntrials, 1);    % store the obtained reward

Qt(1,:)  = Q0;              % initalise Q values

% value simulation
for t = 1:ntrials
    
    PA(t)   = 1./(1+exp(-inv_temp.*(Qt(t,2)-Qt(t,1))));
    ch(t)   = 1 + double(rand()<PA(t));
    
    % outcome
    r(t) = O(t,ch(t));
    
    % compute prediction error
    PE(t) = r(t) - Qt(t,ch(t));
    
    % update value
    Qt(t+1,ch(t)) = Qt(t,ch(t)) + alpha.*PE(t);     % column ch(t) = chosen (1 or 2)
    Qt(t+1,3-ch(t)) = Qt(t,3-ch(t));                % column 3-ch(t) = unchosen (2 or 1)
    
end

%% ------------------------- PLOTS  ---------------------------------------
figure

% display the simulation results for predicted value
subplot(2, 1, 1)
hold on
plot(Qt(:,1), '-b');
plot(Qt(:,2), '-r');
xlim([0, ntrials]);
ylim([0, 1.1]);
xlabel('trials');
ylabel('Q-vvalues');
legend('B','A');

% display the simulation results for prediction error
subplot(2, 1, 2);
hold on
plot(PA, '-k');
plot(ch-1, 'o');
xlim([0, ntrials]);
ylim([-0.1, 1.1]);
xlabel('trials');
ylabel('Proba(choose A)');
legend('proba','choice');
