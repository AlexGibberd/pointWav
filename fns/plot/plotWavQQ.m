function [ qq_theory ] = plotWavQQ( SpecArray,L,coord,varargin )
%PLOTWAVQQ Calculates the appropriate quantities to generate a
% QQplot at any point in rescaled scale time, (a,z)
% Uses a function to look up the nearest estimated point in (a,z)
% Assumes the distribution at this point is Goodman and calculates
% theoretical values for the observed input points
%
%   Inputs:
%       SpecArray - array containing structs 
%           SpecArray(i).Spec.a - sampling points at scale
%           SpecArray(i).Spec.Z - sampling points at Z
%           SpecArray(i).Spec.Coh{a} - estimates of coherence
%       Note: as usual, the number of samples in z at a given a
%       can vary with scale...
%       L - degrees of freedom for goodman
%       coord - scale-time pair (a,z) in (0,1)x(0,1)
%   Optional:
%       ESpec - struct containing the expected spectra
%       Note: this should be structured the same as SpecArray(i),
%           if not given, it defaults at the null model Ecoh=0  
%       dims - [i,j] pair index of coherence to investiage, for two
%              dimensions, we need [1,2] (default)
%       trim - (default 0) sometimes we need to trim large coherences to perform root
%       finding. (see Goodman_QQ_Plots for more detail)
%       wavType - Morlet (default)/Mexican

%% Input parsing
pa = inputParser;  % Need to parse support first
default_ESpec.a = SpecArray(1).a;
default_ESpec.Z = SpecArray(1).Z;
default_ESpec.muC = SpecArray(1).Shat;
default_trim = 0;
default_dims = [1,2];
Na = length(default_ESpec.a);
% Place zeros everywhere... for null model
for j=1:Na
   default_ESpec.muC{j} = zeros(size(SpecArray(1).Shat{j})); 
end
addOptional(pa,'ESpec',default_ESpec,@isstruct);
addOptional(pa,'dims',default_dims,@isnumeric);
addOptional(pa,'trim',default_trim,@isnumeric);
addOptional(pa,'wavType','Morlet');

parse(pa,varargin{:});  
ESpec = pa.Results.ESpec;
dims = pa.Results.dims;
trim = pa.Results.trim;
wavType = pa.Results.wavType;

Nexp = length(SpecArray);

[aidx, zidx] = findClosestSample(ESpec,coord);
obs=[];
for n=1:Nexp
    % Fill observations with closest points in (a,z)
    obs = [obs,SpecArray(n).Coh{aidx}(zidx,dims(1),dims(2))];
    % Do we need to filter for non-compliant coherence here? Nans?
    
end
obs = sort(obs);
% Calculate theoretical points
if(strcmp(wavType,'Morlet'))
    % Complex wavelet
    qq_theory = Goodman_QQ_Plots(ESpec.muC{aidx}(zidx,dims(1),dims(2)),obs,L,trim);
else
    % Real valued wavelet
    qq_theory = Goodman_QQ_Plots_Real(ESpec.muC{aidx}(zidx,dims(1),dims(2)),obs,L,trim);
end
% Plot the result
plot(qq_theory(1:Nexp-trim),obs(1:Nexp-trim),'.');
hold on
plot(linspace(0,1,2),linspace(0,1,2));
title(['QQPlot at position: a=',num2str(ESpec.a(aidx)),', z=',num2str(ESpec.Z{aidx}(zidx))]);
xlabel('theoretical values');
ylabel('observed values');
xlim([0 1]);
ylim([0 1]);

end


