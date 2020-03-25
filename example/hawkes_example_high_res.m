% Script for high-res Hawkes example

load('simulated_hawkes.mat');

%% Pre-processing
Tmin = min(E{1}(1),E{2}(1));
Tmax = max(E{1}(end),E{2}(end));

%% Individual trial analysis
% Wavelet Parameters
wavType = 'Morlet'; % Set the wavelet type for analysis
d = 1;    % Scale/Frequency parameter Morlet. Shape param for Mexican hat
kappa = 10; % Smoothing for twcoh
Na = 100;    % Number of scale levels to sample (100 used in paper)
oversample = 20; % Oversampling for shifts at each scale (20 used in paper)
a_list = linspace(0.05,1,Na);   % List of scales to sample (rescaled between 0 and Amax)

% Precalculate phi to save computation
[phi,lambdas,s,delta] = nystromWav(kappa,'wavType',wavType,'L',10,'Sampling','fixed','NNyst',40);
[dof] = calcDOF(lambdas);

S = tsWP( E, kappa, wavType,'a',a_list,'method','eigenWav',...
    'oversample',oversample,'phi',phi,'lambdas',lambdas,'s',s,'verbose',1 );
     
cu = fzero(@(x) estcoh(x,dof,2,0.95,0 ),[0.000001,0.9999999999]);

plotCWTSpec( S.a,S.Z,S.Coh, 'dims',[1,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 1])

% Fig. 2 - Raster plot of the neuron firing patterns
fh = figure(2);
subplot(2,1,1)
line([E{1},E{1}] , [-0.5 0.5], 'Color', 'k','LineWidth',0.3);
xlabel('Time')
ylabel('$X_1$')
ylim([-0.5,0.5])

subplot(2,1,2)
line([E{2},E{2}], [-0.5 0.5], 'Color', 'k','LineWidth',0.3);
xlabel('Time')
ylabel('$X_2$')
ylim([-0.5,0.5])

