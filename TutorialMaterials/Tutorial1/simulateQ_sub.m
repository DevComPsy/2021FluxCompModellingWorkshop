clear all
close all

%-------------------------------------------------------------------------%
%    SIMULATION OF THE SIMPLE Q-LEARNING MODEL FOR INSTRUMENTAL CONDITIONING
%
%-------------------------------------------------------------------------%

rng('shuffle')
%% -------------------------- INPUTS  -------------------------------------

% NUMBER OF TRIALS
ntrials = 24;
nsess   = 12;
nsub    = 20;


%% ------------------------- PRE_ALLOCATE --------------------------------
Qfull  = nan(ntrials+1, 2,nsub);
PAfull  = nan(ntrials, nsub);
chfull  = nan(ntrials, nsub);

%%
for ksub = 1:nsub
    
    alpha       = rand();
    inv_temp    = 10*rand();
    
    % initialise value and prediction error vectors to store the data
    Qt  = nan(ntrials+1, 2,nsess);
    PE  = nan(ntrials, nsess);
    ch  = nan(ntrials, nsess);
    PA  = nan(ntrials, nsess);
    r   = nan(ntrials, nsess);
    
    for ksess = 1:nsess
        
        %% INITIAL VALUE
        Q0  = [0.5 0.5]; % initial value
        
        %% ------------------------- TASK -----------------------------------
        % Reward history
        O = double([rand(ntrials,1)<0.2, rand(ntrials,1)<0.8]);
        
        % initialise outcome vector
        Qt(1,:,ksess)  = Q0;
        
        % value simulation
        for t = 1:ntrials
            
            PA(t,ksess)   = 1./(1+exp(-inv_temp.*(Qt(t,2,ksess)-Qt(t,1,ksess))));
            ch(t,ksess)   = 1 + (double(rand()<PA(t,ksess)));
            
            % outcome
            r(t,ksess) = O(t,ch(t,ksess));
            
            % compute prediction error
            PE(t,ksess) = r(t,ksess) - Qt(t,ch(t,ksess),ksess);
            
            % update value
            Qt(t+1,ch(t,ksess),ksess) = Qt(t,ch(t,ksess),ksess) + alpha.*PE(t,ksess);
            Qt(t+1,3-ch(t,ksess),ksess) = Qt(t,3-ch(t,ksess),ksess);
            
        end
    end
    
    alphaSim    = alpha;
    betaSim     = inv_temp;
    save(strcat('sim_data_sub',num2str(ksub)),'ch','r','nsess','alphaSim','betaSim')
    
    % for figure
    Qfull(:,:,ksub) = squeeze(mean(Qt,3));
    chfull(:,ksub)  = mean(ch,2);
    PAfull(:,ksub)  = mean(PA,2);
end
% %% ------------------------- PLOTS  ---------------------------------------
figure

% display the simulation results for predicted value
subplot(2, 1, 1)
hold on
errorbar(squeeze(mean(Qfull(:,1,:),3)),squeeze(std(Qfull(:,1,:),0,3))./sqrt(nsub), '-b');
errorbar(squeeze(mean(Qfull(:,2,:),3)),squeeze(std(Qfull(:,2,:),0,3))./sqrt(nsub), '-r');
xlim([0, ntrials]);
ylim([-0.1, 1.1]);
xlabel('trials');
ylabel('Q-vvalues');
legend('B','A');

% display the simulation results for prediction error
subplot(2, 1, 2);
hold on
errorbar(mean(PAfull,2),std(PAfull,0,2)./sqrt(nsub), '-k');
errorbar(mean(chfull-1,2),std((chfull-1),0,2)/sqrt(nsub), 'o');
xlim([0, ntrials]);
ylim([0.4, 1.1]);
xlabel('trials');
ylabel('Proba(choose A)');
legend('proba','choice');
