function [ ] = plotCWTSpec( a,Z,S, varargin )
%PLOTCWTSPEC plots the spectra of the CWT
% Inputs:
%   a = array of sample points in scale
%   Z = cell array of sample points in time at each scale a
%   S = Spectral estimate at each point in (a,Z)
%   (Optional)
%   dims - dimensions to plot in cross-spectra 

% Input parsing
if ndims(S{1})>2
    [~,P,~] = size(S{1});
else
    P=1;
end
Na = length(a);
pa = inputParser;  % Need to parse support first
default_dims = [1,1];

addOptional(pa,'dims',default_dims,@isnumeric);
addOptional(pa,'auxplot',0,@isnumeric);
addOptional(pa,'Amax',1,@isnumeric);
addOptional(pa,'plotFreq',0,@isnumeric);
addOptional(pa,'Delta',1,@isnumeric);
addOptional(pa,'Tmax',1,@isnumeric)

parse(pa,varargin{:});  
dims = pa.Results.dims;
auxplot = pa.Results.auxplot;
Amax = pa.Results.Amax;
plotFreq = pa.Results.plotFreq;
Delta = pa.Results.Delta;
Tmax = pa.Results.Tmax;

% Convert S if multivariate input...
if (P>1)
   for j=1:Na
      S{j} = real(S{j}(:,dims(1),dims(2))); % take first element by default
   end
end

for j=1:Na
    Nz(j)=length(Z{j});
end

t_flat = [];
a_flat = [];
S_flat = [];
for j=1:Na
    t_flat = [t_flat; Z{j}(:)];
    S_flat = [S_flat;squeeze(S{j}(:))];
    for k=1:Nz(j)
       a_flat = [a_flat;a(j)]; 
    end
end
t_flat = t_flat*Tmax;

% Find mean spectra at each scale
for j=1:Na
   mu(j) = nanmean(S{j});
   st(j) = nanstd(S{j});
end
% Interpolate flattened data and plot
[xq,yq] = meshgrid(linspace(0,1,500),linspace(0,Tmax,500));
vq = griddata(a_flat,t_flat,S_flat,xq,yq);

if(auxplot==1)
    subplot(2,2,[1 3]);
end

if(Amax==1)
    scalelabel = 'Normalised Scale a=A/A_T';
else
    scalelabel = 'Scale A';
    a = Amax.*a;    % Change scale for plotting
    a_flat = Amax.*a_flat;
    xq = xq.*Amax;
end

if(plotFreq == 1)
   % If we want to plot in frequencies
   a_flat = Delta./a_flat; 
   xq = Delta./xq;
   scalelabel = 'Frequency (Hz)';
   Amax = max(a_flat);
end

surfc(xq,yq,vq);
shading interp
hold on
plot3(a_flat,t_flat,S_flat,'o')
xlim([0 Amax])
ylim([0 Tmax])
xlabel(scalelabel)
if(Tmax==1)
    ylabel('Normalised Time z=t/T')
else
    ylabel('Time')
end
title('Wavelet Spectra (Real part)');
set(gca,'Ydir','reverse')
view(-90,90) % Set default view from above


if(auxplot==1)
    subplot(2,2,2);
    plot(a,mu);
    xlabel(scalelabel)
    ylabel('Mean')
    subplot(2,2,4);
    plot(a,st);
    xlabel(scalelabel)
    ylabel('Standard deviation')
end
end

