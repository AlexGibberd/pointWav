%% Script to produce the figures found in the accompanying paper 
% To illustrate the method, we consider 6 neuron firing traces from an
% experiment described in "Visual Receptive Field Properties of Neurons 
% in the Mouse Lateral Geniculate Nucleus", Tang et al. 2015 Plos One. The
% data is taken from the supplementary material of that paper.
%
% Specifically we consider cell_id 108 and 117 which indicate high firing
% rates, making them more ameniable to study with pointWav.
load('mouse_subset_sf_108_117.mat');

%% Pre-processing
% First, we pre-process the data so there is a consistent interval to study
Nexp = length(X_1);
Tmax = 0; Tmin = 1000;
for n=1:Nexp
    if(X_1{n}(1) < Tmin || X_2{n}(1) < Tmin )
       Tmin = min(X_1{n}(1),X_2{n}(1));
       shortest = n;
    end
    if(X_1{n}(end) > Tmax || X_2{n}(end) > Tmax)
        Tmax = max(X_1{n}(end),X_2{n}(end));
        longest = n;
    end
end
for n=1:Nexp
    X_1{n} = X_1{n}-Tmin;
    X_2{n} = X_2{n}-Tmin;
end
Tmax = Tmax-Tmin;

%% Individual trial analysis
% As an example, we can consider trial 5 which displays some interesting
% coherent behaviour at larger scales.
E{1} = X_1{5};
E{2} = X_2{5};
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
S = tsWP( E, kappa, wavType,'phi',phi,'lambdas',lambdas,'a',a_list,...
    'method','kernel','oversample',oversample );
% Now eigenvalue
Seig = tsWP( E, kappa, wavType,'a',S.a,'Z',S.Z,'Amax',S.Amax,...
        'method','eigenWav','oversample',oversample,'phi',phi,'lambdas',lambdas,'s',s );
    
% If you want to get 95th percentil of null you can run the below (Asssumes Goodman Coh)
cu = fzero(@(x) estcoh(x,dof,2,0.95,0 ),[0.000001,0.9999999999]);

% Fig. 1 - 3d plot of the estimated wavelet spectra and coherence
% To give comparison between nystrom eigenWavelet approximation and exact
% kernel sampling compare rows: Top) Kernel, Bottom) EigenWavelet (L=10)
figure('Position', [100, 100, 1200, 350])
subplot(2,3,1)
plotCWTSpec( S.a, S.Z, S.Shat, 'dims',[1,1] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
% Note: If you try another example, you may want to change limits here.
colorbar('southoutside')
title('Wavelet Spectra $S_{11}(a,b)$')
subplot(2,3,2)
plotCWTSpec( S.a, S.Z, S.Shat, 'dims',[2,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
colorbar('southoutside')
title('Wavelet Spectra $S_{22}(a,b)$')
subplot(2,3,3)
plotCWTSpec( S.a,S.Z,S.Coh, 'dims',[1,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 1])
colorbar('southoutside')
title('Wavelet Coherence $\gamma_{12}^2(a,b)$')
subplot(2,3,4)
plotCWTSpec( S.a, S.Z, Seig.Shat, 'dims',[1,1] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
% Note: If you try another example, you may want to change limits here.
colorbar('southoutside')
title('Wavelet Spectra $S_{11}(a,b)$')
subplot(2,3,5)
plotCWTSpec( S.a, S.Z, Seig.Shat, 'dims',[2,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
colorbar('southoutside')
title('Wavelet Spectra $S_{22}(a,b)$')
subplot(2,3,6)
plotCWTSpec( S.a,S.Z,Seig.Coh, 'dims',[1,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 1])
colorbar('southoutside')
title('Wavelet Coherence $\gamma_{12}^2(a,b)$')


% Fig. 2 - Raster plot of the neuron firing patterns
fh = figure(2);
for tr=1:Nexp
    subplot(2,1,1)
    if ~isempty(X_1{tr})
        line([X_1{tr};X_1{tr}] , [tr-0.5 tr+0.5], 'Color', 'k','LineWidth',1);
        xlabel('Times (s)')
        ylabel('Trial (cell 108)')
        ylim([0.5,Nexp+0.5])
    end
    subplot(2,1,2)
    if ~isempty(X_2{tr})
        line([X_2{tr};X_2{tr}], [tr-0.5 tr+0.5], 'Color', 'k','LineWidth',1);
        xlabel('Times (s)')
        ylabel('Trial (cell 117)')
        ylim([0.5,Nexp+0.5])
    end
end

%% Multiple Trial Analysis
% If we have multiple trials and we suspect they are drawn i.i.d it is
% reasonable to take an average of the spectra and assess the coherence of
% this.
Smt = [];
for tr=1:Nexp
    E{1} = X_1{tr};
    E{2} = X_2{tr};
    % Estimate spectra (use same sampling points as first example)
    Stem = tsWP( E, kappa, wavType,'a',S.a,'Z',S.Z,'Amax',S.Amax,...
        'verbose',0,'method','kernel','oversample',oversample );
    Smt = [Smt; Stem];  % Store result in indexed struct
end
% To average in pointWav run:
% alpha = percentile to use for empirical confidence interval
Savg = averageSpectra( Smt, 'alpha',0.1 );   

% Fig.3 - plot of the average spectra, and average coherence
figure('Position', [100, 100, 1200, 350])
subplot(1,3,1)
plotCWTSpec( Savg.a, Savg.Z, Savg.mu, 'dims',[1,1] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
colorbar('southoutside')
title('Avg Wavelet Spectra $S_{11}(a,b)$')
subplot(1,3,2)
plotCWTSpec( Savg.a, Savg.Z, Savg.mu, 'dims',[2,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 100])
colorbar('southoutside')
title('Avg Wavelet Spectra $S_{22}(a,b)$')
subplot(1,3,3)
plotCWTSpec( Savg.a, Savg.Z, Savg.muC, 'dims',[1,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 1])
colorbar('southoutside')
title('Avg Wavelet Coherence $\gamma_{12}^2(a,b)$')
% Remark: Although the trial 5 exhibited significant cohernece, this does
% not seem to be the case on average (see reduced coherence in top-right)

% If we assume that trails are i.i.d, then the degress of freedom should be
% 6*dof (the original dof). To access the coherence dervied from Savg.mu,
% call S.Cmu
figure(4)
plotCWTSpec( Savg.a, Savg.Z, Savg.Cmu, 'dims',[1,2] ,'Amax',S.Amax,'Tmax',Tmax);
caxis([0 1])
colorbar('southoutside')
title('Wavelet Coherence (from Avg Spectra) $\gamma_{12}^2(a,b)$')
% The threshold for significance (assuming iid sampling) would then be
% given by
mu_cu = fzero(@(x) estcoh(x,dof*Nexp,2,0.95,0 ),[0.000001,0.9999999999]);


