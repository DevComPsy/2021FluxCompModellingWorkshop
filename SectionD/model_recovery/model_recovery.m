%% Plot model metrics and recoverability %%
% Kate Nussenbaum - katenuss@nyu.edu
% Flux Computational Modeling Workshop
% Last edited: 8/30/21

%clear everything
clear;

%% LOAD DATA %%
load('simulated_data/fits_30-Aug-2021');

%% CREATE CONFUSION MATRICES %%
%-----------------------------------------------------------------%
% STEP 1: Extract best AICs and BICs for each dataset and subject %
%-----------------------------------------------------------------%
%initialize matrices to store AIC and BIC values
best_models_aic = []; 
best_models_bic = [];

%loop through datasets
for dataset = 1:size(model_fits, 1)
    dataset_aics = []; 
    dataset_bics = [];
    
    %loop through models
    for model = 1:size(model_fits, 2)
        dataset_aics = [dataset_aics, model_fits(dataset, model).results.AIC];
        [~, dataset_aic_index] = min(dataset_aics, [], 2);
        dataset_bics = [dataset_bics, model_fits(dataset, model).results.BIC];
        [~, dataset_bic_index] = min(dataset_bics, [], 2);
    end
    best_models_aic = [best_models_aic, dataset_aic_index];
    best_models_bic = [best_models_bic, dataset_bic_index];
end

%-----------------------------------------------------------%
% STEP 2: Create tables with simulated and recovered models %
%-----------------------------------------------------------%

% Get model names
for dataset = 1:size(model_fits, 2)
    model_name{dataset} = model_fits(dataset, 1).sim_model(5:end);
end

%create tables
aic_table = array2table(best_models_aic, 'VariableNames', model_name);
aic_table_stacked = stack(aic_table, 1:size(model_fits, 2),  'NewDataVariableName','Recovered Model',...
          'IndexVariableName','Simulated Model');
bic_table = array2table(best_models_bic, 'VariableNames', model_name);
bic_table_stacked = stack(bic_table, 1:size(model_fits, 2),  'NewDataVariableName','Recovered Model',...
          'IndexVariableName','Simulated Model');      

%------------------------%
% STEP 3: Create heatmap %
%------------------------%

% Plot 
figure;
subplot(1,2,1)
h = heatmap(aic_table_stacked, 'Simulated Model', 'Recovered Model');
h.YDisplayLabels = model_name;
title('AIC');
set(gca,'FontSize',14)
colorbar off
subplot(1,2,2)
h = heatmap(bic_table_stacked, 'Simulated Model', 'Recovered Model');
h.YDisplayLabels = model_name;
title('BIC');
set(gca,'FontSize',14)
colorbar off

%% PLOT MEAN AND MEDIAN AIC AND BIC VALUES %%
%-------------------------------------------------------------------------------%
% STEP 1: Compute mean and median AIC and BIC values for each dataset and model %
%-------------------------------------------------------------------------------%

%initialize matrices
mean_aic = NaN(size(model_fits));
mean_negloglik = NaN(size(model_fits));
num_params = NaN(size(model_fits));
med_aic = NaN(size(model_fits));
mean_bic = NaN(size(model_fits));
med_bic = NaN(size(model_fits));

for dataset = 1:size(model_fits, 1)
    for model = 1:size(model_fits, 2)
        mean_aic(dataset, model) = mean(model_fits(dataset, model).results.AIC);
        mean_negloglik(dataset, model) = mean(model_fits(dataset, model).results.negloglik);
        num_params(dataset, model) = size(model_fits(dataset, model).results.params, 2);
        med_aic(dataset, model) = median(model_fits(dataset, model).results.AIC);
        mean_bic(dataset, model) = mean(model_fits(dataset, model).results.BIC);
        med_bic(dataset, model) = median(model_fits(dataset, model).results.BIC);
    end 
end

%----------------------------------------------%
% STEP 2: % Plot mean and median AICs and BICs %
%----------------------------------------------%

figure;
subplot(2,2,1)
b = bar(mean_aic, 'EdgeColor','black', 'LineWidth', 1);
set(gca, 'xticklabel', model_name);
set(gca,'FontName','Helvetica','FontSize',16);
xlabel('Simulated Model','FontSize',18);
ylabel('Mean AIC', 'FontSize', 18);
leg = legend(model_name, 'Location', 'northeastoutside');
title(leg,'Fit Model');

subplot(2,2,2)
b = bar(med_aic, 'EdgeColor','black', 'LineWidth', 1);
set(gca, 'xticklabel', model_name);
set(gca,'FontName','Helvetica','FontSize',16);
xlabel('Simulated Model','FontSize',18);
ylabel('Median AIC', 'FontSize', 18);
leg = legend(model_name, 'Location', 'northeastoutside');
title(leg,'Fit Model');

subplot(2,2,3)
b = bar(mean_bic, 'EdgeColor','black', 'LineWidth', 1);
set(gca, 'xticklabel',model_name);
set(gca,'FontName','Helvetica','FontSize',16);
xlabel('Simulated Model','FontSize',18);
ylabel('Mean BIC', 'FontSize', 18);
leg = legend(model_name, 'Location', 'northeastoutside');
title(leg,'Fit Model');

subplot(2,2,4)
b = bar(med_bic, 'EdgeColor','black', 'LineWidth', 1);
set(gca, 'xticklabel',model_name);
set(gca,'FontName','Helvetica','FontSize',16);
xlabel('Simulated Model','FontSize',18);
ylabel('Median BIC', 'FontSize', 18);
leg = legend(model_name, 'Location', 'northeastoutside');
title(leg,'Fit Model');


%%
%---------------------------------------------------------------------------------%
% For Flux Workshop: Plot mean AICs, decomposed into two terms, ignore null model %
%---------------------------------------------------------------------------------%

figure;
subplot(1,3,1)
b = bar(2*mean_negloglik(1, 1:3),'FaceColor', [153,51,102]/255, 'EdgeColor','black', 'LineWidth', 1);
set(gca, 'xticklabel', model_name(1:3));
set(gca,'FontName','Helvetica','FontSize',16);
ylim([min(2*mean_negloglik(1, 1:3)) - 5, max(2*mean_negloglik(1, 1:3)) + 10])
xlabel('Model','FontSize',18);
ylabel('-2ln(L)', 'FontSize', 18);

subplot(1,3,2)
b = bar(2*num_params(1, 1:3), 'FaceColor', [211,211,211]/255, 'EdgeColor','black', 'LineWidth', 1);
set(gca, 'xticklabel', model_name(1:3));
set(gca,'FontName','Helvetica','FontSize',16);
ylim([0, 2*max(num_params(1, 1:3)) + 1])
xlabel('Model','FontSize',18);
ylabel('2k', 'FontSize', 18);

subplot(1,3,3)
b = bar([2*mean_negloglik(1, 1:3); 2*num_params(1, 1:3)]', 'stacked' ,'EdgeColor','black', 'LineWidth', 1,'FaceColor','flat');
b(1).CData = [153,51,102]/255;
b(2).CData = [211,211,211]/255;
set(gca, 'xticklabel', model_name(1:3));
set(gca,'FontName','Helvetica','FontSize',16);
ylim([min(2*mean_negloglik(1, 1:3)) - 5, max(2*mean_negloglik(1, 1:3)) + 10])
xlabel('Model','FontSize',18);
ylabel('AIC', 'FontSize', 18);







