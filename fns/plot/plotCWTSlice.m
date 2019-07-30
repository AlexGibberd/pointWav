function [ ] = plotCWTSlice( a,Z,S, varargin )
%PLOTCWTSLICE Plots a slice across scale of cross spectra
% Inputs:
%   a = array of sample points in scale
%   Z = cell array of sample points in time at each scale a
%   S = Spectral estimate at each point in (a,Z)
%   (Optional)
%   Serror - symmetric error bars
%   dims - dimensions to plot in cross-spectra 
%   zplot - position on z (0,1) that we want to take cross specta (default 0.5)
%   lq - lower quartiles
%   uq - upper quartile
%   freq - (default 0) whether to plot against scale or frequency.
%   Amax - default 1. Scales the scale/frequencies

% Input parsing
if ndims(S{1})>2
    [~,P,~] = size(S{1});
else
    P=1;
end
Na = length(a);
pa = inputParser;  % Need to parse support first
default_dims = [1,1];
addOptional(pa,'zplot',0.5,@isnumeric);
addOptional(pa,'lq',[],@iscell);
addOptional(pa,'uq',[],@iscell);
addOptional(pa,'dims',default_dims,@isnumeric);
addOptional(pa,'freq',0,@isnumeric);
addOptional(pa,'Serror',0,@iscell);
addOptional(pa,'Amax',1,@isnumeric);


parse(pa,varargin{:});  
dims = pa.Results.dims;
zplot = pa.Results.zplot;
lq = pa.Results.lq;
uq = pa.Results.uq;
freq = pa.Results.freq;
Serror = pa.Results.Serror;
Amax = pa.Results.Amax;

A = a.*Amax;  % Create absolute scale 
% Convert S if multivariate input...
if (P>1)
   for j=1:Na
      S{j} = real(S{j}(:,dims(1),dims(2))); % take first element by default
   end
end

% For each scale level find sample point closest to zplot

for j=1:Na
    Zj = Z{j};
    [~, idz(j)] = min(abs(Zj-zplot));
    mu(j) = S{j}(idz(j));
    if(iscell(Serror))
        std(j) = Serror{j}(idz(j));
    end
    if (~isempty(lq) && ~isempty(uq))
        lq_plot(j) = mu(j)-lq{j}(idz(j),dims(1),dims(2));
        uq_plot(j) = uq{j}(idz(j),dims(1),dims(2))-mu(j);
    end
end

if(freq == 1)
   A = (1./A).*(2*pi);  % Need to check this is 4pi...
end
if(~isempty(lq) && ~isempty(uq))
    errorbar(A,mu,lq_plot,uq_plot);
elseif(iscell(Serror))
    errorbar(A,mu,std);
else
    plot(A,mu);
end

