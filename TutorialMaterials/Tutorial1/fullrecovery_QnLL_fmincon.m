clear all
close all
clc


%-------------------------------------------------------------------------%
%    SIMULATION OF THE SIMPLE Q-LEARNING MODEL FOR INSTRUMENTAL CONDITIONING
%
%-------------------------------------------------------------------------%

nfl = dir('sim_data_sub*');
nsub = length(nfl);

estimP  = NaN(nsub,2);
simuP   = NaN(nsub,2);

for ksub = 1:nsub
    
    load(strcat('sim_data_sub',num2str(ksub)))
    
    %% -------------------------- INPUTS  -------------------------------------
    
    x0 = [5 0.5];
    xmin = [0 0];
    xmax = [50 1];
    
    options = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000); % These increase the number of iterations to ensure the convergence
    
    [parameters,nll,~,~,~]       = fmincon(@(x) Qlearner(x,ch,r),x0,[],[],[],[],xmin,xmax,[],options);
    
    
    estimP(ksub,:) = parameters;
    simuP(ksub,:) = [betaSim,alphaSim];
    
    
end


figure
subplot(1,2,1)
hold on
plot([0,10],[0,10],'--k')
[b,stats]    = robustfit(simuP(:,1),estimP(:,1));
XX             = linspace(min(simuP(:,1)),max(simuP(:,1)),1000);
[Yf]           = glmval(b,XX,'identity');
XXX            = sortrows([XX',Yf],1);
hModel         =   plot(XXX(:,1),XXX(:,2),'-',...
    'Color',.85*[0,0,1],...
    'LineWidth',2);

plot(simuP(:,1),estimP(:,1),'o',...
    'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0,0,0])
xlabel('true \beta')
ylabel('recovered \beta')

subplot(1,2,2)
hold on
plot([0,1],[0,1],'--k')
[b,stats]    = robustfit(simuP(:,2),estimP(:,2));
XX             = linspace(min(simuP(:,2)),max(simuP(:,2)),1000);
[Yf]           = glmval(b,XX,'identity');
XXX            = sortrows([XX',Yf],1);
hModel         =   plot(XXX(:,1),XXX(:,2),'-',...
    'Color',.85*[0,0,1],...
    'LineWidth',2);
plot(simuP(:,2),estimP(:,2),'o',...
    'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0,0,0])
xlabel('true \alpha')
ylabel('recovered \alpha')

%%
%%==========================================================


function nLL = Qlearner(params,choice,reward)

inv_temp = params(1);
alpha = params(2);

nsess = size(choice,2);
ntrials = size(choice,1);

PA = NaN(ntrials,nsess);
lik = NaN(ntrials,nsess);
Qt = NaN(ntrials,2,nsess);
PE = NaN(ntrials,nsess);

% INITIAL VALUE
Q0  = [0.5 0.5]; % initial value

for ksess = 1:nsess
    
    %% ------------------------- Q SIMULATION --------------------------------
    Qt(1,:,ksess)  = Q0;
    
    % value simulation
    for t = 1:ntrials
        
        PA(t,ksess)   = 1./(1+exp(-inv_temp.*(Qt(t,2,ksess)-Qt(t,1,ksess))));
        
        if choice(t,ksess) == 1
            lik(t,ksess) = 1 - PA(t,ksess);
        elseif choice(t,ksess) == 2
            lik(t,ksess) = PA(t,ksess);
        end
        
        % compute prediction error
        PE(t,ksess) = reward(t,ksess) - Qt(t,choice(t,ksess),ksess);
        
        % update value
        Qt(t+1,choice(t,ksess),ksess) = Qt(t,choice(t,ksess),ksess) + alpha.*PE(t,ksess);
        Qt(t+1,3-choice(t,ksess),ksess) = Qt(t,3-choice(t,ksess),ksess);
        
    end
end
nLL = -sum(log(lik(:)));

end