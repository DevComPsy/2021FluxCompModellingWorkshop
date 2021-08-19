%% Plot model recoverability

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


%% Make heatmaps
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
