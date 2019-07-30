function [Q] = plotThresh(S,dof,varargin)
%PLOTTHRESH Plots proportion of samples above a coherence threshold relative to
% percentile of Goodman cdf

%   Inputs:
%       SpecArray - array containing structs 
%           S(i).a - sampling points at scale
%           S(i).Z - sampling points at Z
%           S(i).Coh{a} - estimates of coherence
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
%       wavType - Morlet (default)/Mexican
%       alpha - percentile at which to threshold

%% Input parsing
pa = inputParser;  % Need to parse support first
default_ESpec.a = S(1).a;
default_ESpec.Z = S(1).Z;
default_ESpec.muC = S(1).Shat;
default_dims = [1,2];
default_alpha = 0.95;
Na = length(default_ESpec.a);
% Place zeros everywhere... for null model
for j=1:Na
   default_ESpec.muC{j} = zeros(size(S(1).Shat{j})); 
end
addOptional(pa,'ESpec',default_ESpec,@isstruct);
addOptional(pa,'dims',default_dims,@isnumeric);
addOptional(pa,'wavType','Morlet');
addOptional(pa,'alpha',default_alpha,@isnumeric);

parse(pa,varargin{:});  
ESpec = pa.Results.ESpec;
dims = pa.Results.dims;
alpha = pa.Results.alpha;
wavType = pa.Results.wavType;

Nexp = length(S);
Q = cell([1,Na]);
[~,P,~] = size(cell2mat(S(1).Shat(1)));

for j=1:Na
    % Get number of shifts at this scale
    Kj = length(cell2mat(S(1).Z(j)));
    Q{j} = zeros([Kj,1]);
    for k=1:Kj
        Expected = ESpec.muC{j}(k,dims(1),dims(2));
        cu = fzero(@(x) estcoh(x,dof,P,alpha,Expected ),[0.000001,0.9999999999]);
        p = 0; % index for counting items above threshold
        for n=1:Nexp
            Chats = cell2mat(S(n).Coh(j));
            Chat = Chats(k,dims(1),dims(2));
            if(Chat>cu)
                p = p+1;
            end
        end
        Q{j}(k) = p/Nexp;
    end
end

plotCWTSpec( S(1).a,S(1).Z,Q)

end

