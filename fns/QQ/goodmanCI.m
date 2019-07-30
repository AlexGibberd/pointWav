function [cl,cu,gam] = goodmanCI(alphal,alphau,N,P,L)
%% Generates approximate confidence interval over unit square for coherence

gam = linspace(0,0.9999,N); % Investigate grid of true params
for i=1:length(gam)
    cl(i)= fzero(@(x) estcoh(x,L,P,alphal,gam(i)),[0.000001,0.9999999999]);
    cu(i)= fzero(@(x) estcoh(x,L,P,alphau,gam(i)),[0.000001,0.9999999999]);
end

end