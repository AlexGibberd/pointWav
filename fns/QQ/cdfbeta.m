function h = cdfbeta(n,p,x,gam)
%CDFBETA replaces cdfcohfin if 
p = 0;
gam =0; % Only works for true coherence being zero..

h = betacdf(x,0.5,n/2 - 1);

end

