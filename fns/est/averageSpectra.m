function [ Savg ] = averageSpectra( Spec, varargin )
%AVERAGESPECTRA computes empirical average of spectra (given they are
% the same size)
% INPUTS:
%   Spec - Nexp array of structs Spec(n)
%       Spec(n).a
%       Spec(n).Amax
%       Spec(n).Z
%       Spec(n).Shat
%       Spec(n).Coh (Optional)
% (Optional)
%   alpha - level to calculate lower and upper quanitles at.. (defauly 0.05)
%   Strue - strucutre array with true spectra for all relevent scale/shift
%   samples (default empty)
%
% Alex Gibberd 30/05/2019


pa = inputParser;  % Need to parse support first
default_alpha = 0.05;
default_Strue = [];
addOptional(pa,'alpha',default_alpha,@isnumeric);
addOptional(pa,'Strue',default_Strue,@isstruct);

parse(pa,varargin{:});  
alpha = pa.Results.alpha;
Strue = pa.Results.Strue;

Nexp = length(Spec);
% Replicate sample point
Savg.Z = Spec(1).Z;
Savg.a = Spec(1).a;
Savg.Amax = Spec(1).Amax;
Na = length(Spec(1).a);

% Check equal dimensionality for all experiments
if (ndims(Spec(1).Shat{1})>2)
    %P>1 case
    for j=1:Na
        [T(j),P,~] = size(Spec(1).Shat{j});
    end
else
    % P=1 case (univariate)
    P=1;
    for j=1:Na
        T(j) = length(Spec(1).Shat{j});
    end
end

if (P>1)
    for n=1:Nexp
        for j=1:Na
            [K,~,~] = size(Spec(n).Shat{j});    % Get number of timepoints at this scale
            % On first experiment init arrays
            if(n==1)
                Sl{j} = zeros([Nexp,K,P,P]);
                Cohl{j} = zeros([Nexp,K,P,P]);
            end
            for k=1:K
               Sl{j}(n,k,:,:) = squeeze(Spec(n).Shat{j}(k,:,:));
               Cohl{j}(n,k,:,:) = squeeze(Spec(n).Coh{j}(k,:,:));
            end
        end
    end
    for j=1:Na
       Sav{j} = squeeze(nanmean(Sl{j},1));    % compute mean
       Sstd{j} = squeeze(nanstd(Sl{j},0,1));   % Compute std
       Cav{j} = squeeze(nanmean(Cohl{j},1));    % compute mean
       Cstd{j} = squeeze(nanstd(Cohl{j},0,1));   % Compute std
       Slq{j} = squeeze(quantile(Sl{j},alpha,1));    % compute lower quanitle
       Suq{j} = squeeze(quantile(Sl{j},1-alpha,1));   % compute upper quantile
       Clq{j} = squeeze(quantile(Cohl{j},alpha,1));    % compute lower quanitle
       Cuq{j} = squeeze(quantile(Cohl{j},1-alpha,1));   % compute upper quanitle
       Ctem = zeros(size(Sav{j})); 
       for k=1:size(Sav{j},1)
            for l=1:P
                for m=1:P
                    % Should be real..
                    Ctem(k,l,m) = abs(Sav{j}(k,l,m))^2./(Sav{j}(k,l,l)*Sav{j}(k,m,m)) ;
                end
            end
        end
        Cmu{j} = Ctem;

    end
    Savg.mu = Sav;
    Savg.std = Sstd;
    Savg.lq = Slq;
    Savg.uq = Suq;
    Savg.muC = Cav;
    Savg.stdC = Cstd;
    Savg.lqC = Clq;
    Savg.uqC = Cuq;
    Savg.Cmu = Cmu;
else
    %% Univariate case...
    for n=1:Nexp
        for j=1:Na
            if(n==1)
                Sl{j} = zeros([Nexp,size(Spec(1).Shat{j},1)]);
            end
            Sl{j}(n,:) = Spec(n).Shat{j};
        end
    end
    for j=1:Na
       Sav{j} = mean(Sl{j}(:,:),1);    % compute mean
       Sstd{j} = std(Sl{j}(:,:),0,1);   % Compute std
       Slq{j} = quantile(Sl{j}(:,:),alpha,1);    % compute lower quantile
       Suq{j} = quantile(Sl{j}(:,:),1-alpha,1);   % Compute upper quantile
    end
    Savg.mu = Sav;
    Savg.std = Sstd;
    Savg.lq = Slq;
    Savg.uq = Suq; 
end

% If Strue is given then calculate error
if(isstruct(Strue))
    for j=1:Na
        Savg.err{j} = abs(Savg.mu{j}-Strue.mu{j}');
    end
end




end

