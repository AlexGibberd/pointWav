%% Script to produce the figures found in the accompanying paper 
% This example uses pre-simulated data from a bi-variate Hawkes process

load('simulated_hawkes.mat');

%% Pre-processing
Tmin = min(E{1}(1),E{2}(1));
Tmax = max(E{1}(end),E{2}(end));

%% Individual trial analysis
% Wavelet Parameters
wavType = 'Morlet'; % Set the wavelet type for analysis
d = 1;    % Scale/Frequency parameter Morlet. Shape param for Mexican hat
kappa = 10; % Smoothing for twcoh
Na = 20;    % Number of scale levels to sample (100 used in paper)
oversample = 10; % Oversampling for shifts at each scale (20 used in paper)
a_list = linspace(0.05,1,Na);   % List of scales to sample (rescaled between 0 and Amax)

% Precalculate phi to save computation
[phi,lambdas,s,delta] = nystromWav(kappa,'wavType',wavType,'L',10,'Sampling','fixed','NNyst',40);
[dof] = calcDOF(lambdas);

% Calculate the spectra first kernel sampling
S = tsWP( E, kappa, wavType ,'a',a_list,...
    'method','kernel','oversample',oversample,'verbose',1 );
% Now eigenvalue
% Seig = tsWP( E, kappa, wavType,'a',S.a,'Z',S.Z,'Amax',S.Amax,...
%         'method','eigenWav','oversample',oversample,'phi',phi,'lambdas',lambdas,'s',s );
    
% If you want to get 95th percentil of null you can run the below (Asssumes Goodman Coh)
cu = fzero(@(x) estcoh(x,dof,2,0.95,0 ),[0.000001,0.9999999999]);

% Fig. 1 - 3d plot of the estimated wavelet spectra and coherence
% To give comparison between nystrom eigenWavelet approximation and exact
% kernel sampling compare rows: Top) Kernel, Bottom) EigenWavelet (L=10)
figure('Position', [100, 100, 1200, 350])
subplot(1,3,1)
plotCWTSpec( S.a, S.Z, S.Shat, 'dims',[1,1] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
% Note: If you try another example, you may want to change limits here.
colorbar('southoutside')
title('Wavelet Spectra $S_{11}(a,b)$')
subplot(1,3,2)
plotCWTSpec( S.a, S.Z, S.Shat, 'dims',[2,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
colorbar('southoutside')
title('Wavelet Spectra $S_{22}(a,b)$')
subplot(1,3,3)
plotCWTSpec( S.a,S.Z,S.Coh, 'dims',[1,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 1])
colorbar('southoutside')
title('Wavelet Coherence $\gamma_{12}^2(a,b)$')



