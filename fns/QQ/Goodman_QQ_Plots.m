function theoretical_values = Goodman_QQ_Plots(mean,observed_values,kk,trim)

L = length(observed_values);

pp=(1:L)./(L+1);

% find quantile of Goodman dist for coherence
% for the AR rho^2=0.36 for all frequencies!

theoretical_values = zeros(1,L);

% Why do we sometimes need to not go through the entire list here. Results
% in error "The function values at the interval endpoints must differ in sign."
% Removing some entries, i.e. L-trim seems to enable the procedure to run.
for j=1:L-trim
    prob=pp(j);
    theoretical_values(j) = fzero(@(x) estcoh(x,kk,2,prob,mean),[0.000001,0.999]);
end