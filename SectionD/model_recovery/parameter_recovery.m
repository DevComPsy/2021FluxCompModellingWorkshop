%% Plot parameter recoverability %%
% Kate Nussenbaum - katenuss@nyu.edu
% Flux Computational Modeling Workshop
% Last edited: 8/30/21

clear;

%% LOAD SIMULATED DATA AND MODEL FITS %%

%simulated data
load('simulated_data/data_30-Aug-2021');

%model fits
load('simulated_data/fits_30-Aug-2021');

%% GET SIMULATED AND RECOVERED PARAMETERS FROM EACH MODEL %%
for m = 1:length(sim_data)
    sim_params = zeros(size(model_fits(m,m).results.params));
    fit_params = zeros(size(model_fits(m,m).results.params));
    for s = 1:length(sim_data(1).sub_data)
        if ~isempty(sim_data(m).sub_data(s).params)
            sim_params(s, :) = sim_data(m).sub_data(s).params;
            fit_params(s, :) = model_fits(m, m).results.params(s, :);
        end
    end
    params(m).sim = sim_params;
    params(m).fit = fit_params;
end

%% PLOT CORRELATIONS BETWEEN SIMULATED AND RECOVERED PARAMETERS %%

for m = 1:length(params)
    if ~isempty(params)
        figure; % new figure for each model
        n_params = sim_data(m).n_params; %get number of parameters
        param_names = sim_data(m).param_names; %get parameter names
        
        for p = 1:n_params
            subplot(1, n_params, p); %new subplot for each parameter
            
            %make scatter plot
            scatter(params(m).sim(:, p), params(m).fit(:, p), 10, 'MarkerFaceColor', [102 102 255]./255, 'MarkerEdgeColor', 'k');
            
            %add linear regression line
            hold on;
            P = polyfit(params(m).sim(:, p), params(m).fit(:, p), 1);
            yfit = P(1)*params(m).sim(:, p)+P(2);
            plot(params(m).sim(:, p),yfit,'k-.', 'LineWidth', 3);
            
            %add parameter name to plot
            title([param_names(p), round(P(1), 2)], 'Interpreter', 'none');
            set(gca, 'FontSize', 12);
        end
        
        %add model name to plot
        sgtitle(sim_data(m).function(5:end), 'FontSize', 20);
    end
    
end

