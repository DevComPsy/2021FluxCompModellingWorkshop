function val = log_lik(theta,cdf,C1,C2)
%LOG_LIK Calculate negative log likelihood
%   Inputs: theta, cdf, weights
g=cdf(theta);
val=-sum(C1.*log(g)+C2.*log(1-g));
end

