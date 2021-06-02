% function [policy] = softmax(Vs,params)
% 
% Softmax function to transform choice values into choice policy
% 
% input:
% @Vs: array of choice values for all available options at trial t
% @params.tau: decision temperature parameter (cave: the larger, the more random - inverse of the 'inverse temperature')
% 
% output:
% @policy: array of choice probabilities for each available option
% 
% Tobias Hauser, 06/2021
%  
function [policy] = softmax(Vs,params)

%% get parameter
tau = params.tau;


%% calculate policy

% remove max to avoid numerical overflow (cf Friston spm_softmax)
% already divide by temperature tau
VT = (Vs - max(Vs)) / tau; 

% softmax transformation 
policy = (exp(VT) / sum(exp(VT)));