function [dof] = calcDOF(lambdas)
%CALCDOF Calculates DOF given a set of lambdas

% Assume that all etas are invluded
L = length(lambdas);
etas = (L/norm(1./lambdas,1)).*(1./lambdas);

dof = (sum(etas).^2)/sum(etas.^2);
    
end

