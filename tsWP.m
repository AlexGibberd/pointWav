function [ Spec ] = tsWP( E, kappa, wavType, varargin )
% tsWP - Calculates the temporally smoothed wavelet periodogram
% for a p dimensional event stream. Evaluates the spectra at a set of
% points in the viable time-scale triangle.
%
% INPUTS:
%   E - cells of event stream (possibly multivariate) E{1}, E{2}, ..., E{p}
%   kappa - smoothing level relative to A
%   wavType - type of wavelet (default:'Morlet')
%
% OPTIONAL:
%   method - type of method to use, either "eigenWav" (default), or "kernel"
%   phi - Nystrom approximation of eigenwavelets *If precomputed*
%   lambdas - eigenvalues (1/etas) for approximation *If precomputed*
%   s - Nystrom approximation points *If precomputed*
%   NNyst - Number of nystrom samples to use (if creating approx)
%   L - Number of eigenwavelets to use
%   d - frequency/scale link parameter for Morlet (default d=1)
%   a - array of scales to sample a = [1:Na] in (0,1)
%   Amax - upper limit for scale levels (absolute) default T/2*SupScale
%   Z - cell of time sample positions Z{Na}[1:Nz(a)]
%   verbose - 1/0 (default 0)
%
% OUTPUT:
%   S - estiamte of univariate wavelet spectra
%       This is stored as a struct
%       S.a = a, S.Z = Z, 
%       S.Shat[idx_a,idx_t] = estimated spectra
%
%   Alex Gibberd 30/07/2019

    % Input parsing
    pa = inputParser;  % Need to parse support first

    % Is input uni-or multivariate
    if(iscell(E)==1)
        [~,P] = size(E);
        T = findMaxT(E);
    else
        P=1;
        T = E(end); % Set time period to be last event
    end

    % Work in relative scales a in (0,1] 
    a_l = 0.05;  % Lower bound for a..can we get a theoretically motivated number here?
    Na = 10;   % Number of scale scamples
    default_a = linspace(a_l,1,Na);
    alpha = 1e-10;  % For calculating approx support.
    supK = approxSup(wavType,alpha,kappa);  % Support of eigenwavelets at A=1
    default_Amax = T/(2*supK);

    addOptional(pa,'d',1,@isnumeric);
    addOptional(pa,'a',default_a,@isnumeric);
    addOptional(pa,'Amax',default_Amax,@isnumeric);
    addOptional(pa,'NNyst',200,@isnumeric);
    addOptional(pa,'phi',0,@isnumeric);
    addOptional(pa,'lambdas',0,@isnumeric);
    addOptional(pa,'s',0,@isnumeric);
    addOptional(pa,'Z',cell(1),@iscell);
    addOptional(pa,'L',10,@isnumeric);
    addOptional(pa,'verbose',0,@isnumeric);
    addOptional(pa,'method','eigenWav');
    addOptional(pa,'oversample',3);


    % Parse optional inputs
    parse(pa,varargin{:});  

    Z = pa.Results.Z;
    d = pa.Results.d;
    a = pa.Results.a;
    Amax = pa.Results.Amax;
    NNyst = pa.Results.NNyst;
    phi = pa.Results.phi;
    lambdas = pa.Results.lambdas;
    s = pa.Results.s;
    L = pa.Results.L;
    verbose = pa.Results.verbose;
    method = pa.Results.method;
    oversample = pa.Results.oversample;
    
    Na = length(a); % Get number of scales

    % If no approximate wavelet is passed then create
    if(strcmp(method,'eigenWav'))
        if (length(phi) == 1 || length(lambdas) == 1 || length(s) == 1)
            disp('Creating Nystrom Eigenwavelets');
            [phi,lambdas,s,delta] = nystromWav(kappa,'wavType',wavType,...
                'L',L,'Sampling','fixed','alpha',alpha,'NNyst',NNyst);
        end
    end

    % Construct Z if required
    if(isempty(Z{1})==1)
        Z = constructZ(a,Amax,supK,T,oversample);
    end
    % Get the number of samples at each scale
    for j=1:Na
        Nz(j) = length(Z{j});
    end 

    % Set up output structure
    Spec = struct; 
    Spec.a = a;
    Spec.Amax = Amax;
    Spec.Z = Z;

    if (P==1)
        Spec.Shat = cell([Na,1]);
    else
        % Store cross spectra and coherehnce
        Spec.Shat = cell([Na,1]);
        Spec.Coh = cell([Na,1]);
    end

    % For each scale and central wavelet position
    for j=1:Na
        if(verbose == 1)
            disp(['Calculating scale: ',int2str(j),'/',int2str(Na)]);
        end
        % Init the solution array..
        Spec.Shat{j} = squeeze(nan([Nz(j),P,P]));
        Spec.Coh{j} = squeeze(nan([Nz(j),P,P]));
        for k=1:Nz(j)
            % Min and max absolute position of support
            b0 = Z{j}(k);
            %   Since we have already created phi we don't need to 
            %   bother about wavType here..
            if(strcmp(method,'eigenWav'))
               [S,C] = wavShat(E,a(j)*Amax,b0,phi,lambdas,s,kappa,wavType,supK,L);
            else
               [S,C] = kernelShat(E,a(j)*Amax,b0,kappa,supK,wavType); 
            end
                % Package the results
            if(P==1)
               Spec.Shat{j}(k) = real(S);  % Just to make sure take real
            else
               Spec.Shat{j}(k,:,:) = S;
               Spec.Coh{j}(k,:,:) = C;
            end
        end
    end

end

%% AUXILARY FUNCTIONS BELOW...

function T = findMaxT(E)
     % Find last event across all streams
    [~,P] = size(E);
    T=0;
    for j=1:P
        T = max(E{j}(end),T);
    end
end


function [ S,C ] = kernelShat(E,A,b,kappa,supK,wavType)
    
    [S,C]  = SKernel(E,A,b,wavType,kappa,supK);

end

function [ S,C ] = wavShat(E,A,b,phi,lambdas,s,kappa,wavType,supK,L)
%wavShat - a function to actually estimate S at a given scale/shift

    % E isnt always a cell. Only is cell in case of mv PP
    if(iscell(E))
        P = length(E);
        Ti = [];
        for i=1:P
            Ti = [Ti;E{i}(end)];
        end
        T = max(Ti);
    else
        % In 1d operate on vectors..
        T = E(end);
        P = 1;
    end
    B = b*T;    % Absolute shift
    
    % Scale the eigenvalues
    etas = (L/norm(1./lambdas,1)).*(1./lambdas);

    if(P>1)
        % Multivariate 
        for i=1:P
            % We get eigenwav transform (note option)
            vj(i,:) = eigenWav(phi,lambdas,s,E{i},kappa,A,B,'wavType',...
            wavType,'SupK',supK,'padding',0,'transform',1);
        end
        S = zeros([P,P]);   % Init
        for l=1:L
            S = S + etas(l)*(vj(:,l)*ctranspose(vj(:,l)));
        end
        % We then average across tapers
        S = S./L;
        for i=1:P
            for j=1:P
                % Should be real..
                C(i,j) = abs(S(i,j))^2./(S(i,i)*S(j,j)) ;
            end
        end

    else
        % Univariate
        vj = eigenWav(phi,lambdas,s,E,kappa,A,B,'wavType',...
            wavType,'SupK',supK,'transform',1);
        S = 0;
        for l=1:L
            S = S + etas(l)*vj(l)*conj(vj(l));
        end
        S = S./L;
        C = 0;  % No coherence in 1d case
    end
    
end



