%% Plot model recoverability

%Get data
ntrl = cfg_sim.ntrl;
nsubs = cfg_sim.nsims;

% First, compute AIC for null model
null_nll = -1*sum(log(.5*ones(ntrl, 1)));
null_AICs = 2*null_nll*ones(nsubs, 1);

% Concatenate AICs for each dataset
oneLR_AICs = [out_1LRdata_1LR(:, 2), out_1LRdata_2LR(:, 2), null_AICs];
twoLR_AICs = [out_2LRdata_1LR(:, 2), out_2LRdata_2LR(:, 2), null_AICs];
nullModel_AICs = [out_nulldata_1LR(:, 2), out_nulldata_2LR(:, 2), null_AICs];

%Find best fitting model
[oneLR_bestAIC,  oneLR_best] = min(oneLR_AICs, [], 2);
[twoLR_bestAIC,  twoLR_best] = min(twoLR_AICs, [], 2);
[nullModel_bestAIC,  nullModel_best] = min(nullModel_AICs, [], 2);

%Make table
sim_model = [ones(nsubs, 1); 2*ones(nsubs,1); 3*ones(nsubs, 1)];
fit_model = [oneLR_best; twoLR_best; nullModel_best];
model_table = table(sim_model, fit_model);

%% Find best-fitting model for 1 LR data

%Compute mean and median AICs
mean_AICs = mean(oneLR_AICs);
med_AICs = median(oneLR_AICs);
se_AICs = std(oneLR_AICs)/sqrt(nsubs);

%Find proportion of subs best fit by each model
prop_best_fit = [sum(oneLR_best == 1), sum(oneLR_best == 2), sum(oneLR_best == 3)]./nsubs;

figure;
subplot(1,3,1);
b = bar(mean_AICs, 'FaceColor', [153,51,102]/255, 'EdgeColor', 'black', 'LineWidth', 1);
set(gca, 'xticklabel',{'1 LR', '2 LR', 'Null'});
set(gca,'FontName','Helvetica','FontSize',16);
ylim([110 140])
xlabel('Model','FontSize',18);
ylabel('Mean AIC', 'FontSize', 18);

subplot(1,3,2);
b = bar(med_AICs, 'FaceColor', [153,51,102]/255, 'EdgeColor','black', 'LineWidth', 1);
set(gca, 'xticklabel',{'1 LR', '2 LR', 'Null'});
set(gca,'FontName','Helvetica','FontSize',16);
ylim([110 140])
xlabel('Model','FontSize',18);
ylabel('Median AIC', 'FontSize', 18);

subplot(1,3,3);
b = bar(prop_best_fit, 'FaceColor', [153,51,102]/255, 'EdgeColor','black', 'LineWidth', 1);;
set(gca, 'xticklabel',{'1 LR', '2 LR', 'Null'});
set(gca,'FontName','Helvetica','FontSize',16);
ylim([0 1])
xlabel('Model','FontSize',18);
ylabel('Proportion of Participants Best Fit', 'FontSize', 18);


%% Make heatmap showing model recoverability
figure;
h = heatmap(model_table, 'sim_model', 'fit_model');
h.Colormap = winter;
h.Title = 'Model Recoverability';
h.XLabel = 'Simulated Model';
h.YLabel = 'Recovered Model';
h.XDisplayLabels = ["1 LR" "2 LR" "Null"];
h.YDisplayLabels = ["1 LR" "2 LR" "Null"];
set(gca,'FontName','Helvetica','FontSize',16);

saveas(gcf, heatmap_name);
