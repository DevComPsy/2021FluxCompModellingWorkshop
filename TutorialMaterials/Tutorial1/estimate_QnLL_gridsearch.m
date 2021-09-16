clear all
close all
clc


%-------------------------------------------------------------------------%
%    SIMULATION OF THE SIMPLE Q-LEARNING MODEL FOR INSTRUMENTAL CONDITIONING
%
%-------------------------------------------------------------------------%


load('sim_data_sess')

%% -------------------------- INPUTS  -------------------------------------
% NUMBER OF TRIALS
ntrials = 24;

% INITIAL VALUE
Q0  = [0.5 0.5]; % initial value

% PARAMETERS
alpha_mat   = linspace(0,1,201);
temp_mat    = linspace(0,5,201);

nLL = NaN(numel(alpha_mat),numel(temp_mat));

for a = 1:length(alpha_mat)
    alpha = alpha_mat(a);
    
    for b = 1:length(temp_mat)
        inv_temp = temp_mat(b);
        
        
        PA = NaN(ntrials,nsess);
        lik = NaN(ntrials,nsess);
        Qt = NaN(ntrials,2,nsess);
        PE = NaN(ntrials,nsess);
        
        for ksess = 1:nsess
            
            %% ------------------------- Q SIMULATION --------------------------------
            Qt(1,:,ksess)  = Q0;
            
            % value simulation
            for t = 1:ntrials
                
                PA(t,ksess)   = 1./(1+exp(-inv_temp.*(Qt(t,2,ksess)-Qt(t,1,ksess))));
                
                if ch(t,ksess) == 1
                    lik(t,ksess) = 1 - PA(t,ksess);
                elseif ch(t,ksess) == 2
                    lik(t,ksess) = PA(t,ksess);
                end
                
                % compute prediction error
                PE(t,ksess) = r(t,ksess) - Qt(t,ch(t,ksess),ksess);
                
                % update value
                Qt(t+1,ch(t,ksess),ksess) = Qt(t,ch(t,ksess),ksess) + alpha.*PE(t,ksess);
                Qt(t+1,3-ch(t,ksess),ksess) = Qt(t,3-ch(t,ksess),ksess);
                
            end
        end
        
        nLL(a,b) = -sum(log(lik(:)));
    end
end


%%
figure

imagesc(flipud(nLL))
hold on
xlabel('inverse temperature')
ylabel('learning rate')

xt = linspace(1,length(temp_mat),11);
xtl = linspace(min(temp_mat),max(temp_mat),11);

yt = linspace(1,length(alpha_mat),11);
ytl = linspace(min(alpha_mat),max(alpha_mat),11);

set(gca,'XLim',[1 length(temp_mat)],...
    'XTick',xt,...
    'XTickLabel',xtl,...
    'YLim',[1 length(alpha_mat)],...
    'YTick',yt,...
    'YTickLabel',fliplr(ytl))


colorbar


[I,J] = find(nLL == min(min(nLL)));

disp(strcat('alpha = ',num2str(alpha_mat(I))))
disp(strcat('inv. temp = ',num2str(temp_mat(J))))


plot(J,length(alpha_mat)-I,'o',...
    'MarkerFaceColor',[1,0,0],...
    'MarkerEdgeColor',[1,0,0])

text(20,35,strcat('true \alpha = ',num2str(alphaSim)),'Fontsize',14,'color',[1, 1, 1])
text(20,25,strcat('true \beta = ',num2str(betaSim)),'Fontsize',14,'color',[1, 1, 1])


