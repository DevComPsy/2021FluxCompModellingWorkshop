% function [V,PE] = learn_RescWagn(V,outcome,params)
%
% Simple Rescorla-Wagner learning function that updates value V based on a
% reward prediction error (PE) and a learning rate alpha.
%
% inputs:
% @V: Value at trial t
% @outcome: reward/outcome at trial t
% @params: structure containig the learning rate parameter.
%           if this is a single learning rate, field should be called
%           params.alpha
%           for two learning rates, please name parameters as follows
%           params.alpha_pos: positive learning rate
%           params.alpha_ned: negative learning rate
%
% outputs:
% @V: expected value at t+1
% @PE: prediction error at t
%
% Tobias Hauser, 06/2021
%
function [V,PE] = learn_RescWagn(V,outcome,params)


%% calculate prediction error
PE = outcome - V;

%% get learning rate
if PE >= 0 && isfield(params,'alpha_pos')   % separate learning rates
    alpha = params.alpha_pos;
elseif PE < 0 && isfield(params,'alpha_neg')    % separate learning rates
    alpha = params.alpha_neg;
else                                            % only one learning rate
    alpha = params.alpha;
end

%% update value
V = V + alpha * PE;