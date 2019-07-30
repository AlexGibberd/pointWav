function [phi,lambdas,s,delta] = nystromWav(kappa,varargin)
%NYSTROMWAV Creates Nystrom Approximated eigen-wavelets
% Note: All parameterisations are in terms of the mother wavelet (a=1)
% Required:
%   kappa - smoothing window size
% Optional:
%   wavType - 'Morlet' or 'Mexican'
%   L  - number of eigenfunctions to approximate
%   delta - sampling interval (default => N=1000)
%   Sampling - 'fixed' (default), 'normal', 'uniform', the type of sampling for nystrom and quadrature 
%   Seed - seed for random number generation. default 123
%   alpha - proportion of wavelet energy to maintain in support default 1e-7
%   NNyst - Number of nustrom samples

%
% Outputs:
%   phi - matrix (NxL) representing sampled eigen-wavelts 
%   lambdas - vector of eigenvalues 1/etas
%   s - sample (Nystrom) points for approximation
%   delta - spacing for nystrom/quadrature sampling

% First parse the optional inputs
% we decide on the support of the wavelet by looking at 
% cumulative power of the kernel

    default_alpha = 1e-7; % 1- proportion of wavelet energy to maintain
    default_L = 10;
    default_seed = 123;
    pa = inputParser;  % Need to parse support first
    addOptional(pa,'delta',0,@isnumeric);
    addOptional(pa,'wavType','Morlet');
    addOptional(pa,'L',default_L,@isnumeric);
    addOptional(pa,'Sampling','fixed')
    addOptional(pa,'Seed',default_seed,@isnumeric);
    addOptional(pa,'alpha',default_alpha,@isnumeric);
    addOptional(pa,'NNyst',200,@isnumeric);



% Parse first lot of parameters..
    parse(pa,varargin{:});  
    delta = pa.Results.delta;
    wavType = pa.Results.wavType;
    Sampling = pa.Results.Sampling;
    L = pa.Results.L;
    alpha = pa.Results.alpha;
    N = pa.Results.NNyst;

% Seed random number generator
    rng(pa.Results.Seed)

    
% If no delta is set, then estimate delta from l2 norm of K
    if(delta==0)
        SupK = approxSup(wavType,alpha,kappa);
        if(SupK==0)
           disp('Error: wavType not supported (Kernel calc)');
           return 
        end
        delta = SupK/(N/2);
    end
    
% Construct the grid and perform eigen-decomposition
    t = -SupK:delta:SupK;
    Nt = length(t);
    if(strcmp(Sampling,'fixed'))
        % Deterministic sampling
        s = t; % Being explicit here, set Nystrom = Quadrature
    elseif(strcmp(Sampling,'normal'))
        % Normal sampling over real line.
        % Uses SupK to specify Gaussian width
        s = randn([1,Nt]);
        s = sort(s);
%         s = ((2*SupK)./(max(s)-min(s))).*s;
        s = (SupK/1.97).*s;
%         disp('Normal sampling not implemented yet');
%         return;
    elseif(strcmp(Sampling,'uniform'))
        % Uniform sampling
        s = -SupK+(2*SupK).*rand([1,Nt]);
        s = sort(s);
    else
       disp('Unknown sampling type')
       return; 
    end
    [S,T] = meshgrid(s,t);
    if(strcmp(wavType,'Morlet')) 
        % We assume d=1..
        K = errFnKernelmorlet(S,T,kappa);
    elseif(strcmp(wavType,'Mexican'))
        % We assume scale = 1, sigma = 1
        K = mexicanKernel(S,T,1,1,kappa);
    end
    
% Check that there are not too many tapers
    if (L>size(K,1)-1)
        disp('Error: #nystrom samples too small to approximate this many eigenwavelets')
        return;
    end
    
% Perform actual eigen-decomposition
% Note: we assume that w = delta (regular sampling)
    [phi,etas] = eigs(delta*K,L);
    etas = diag(etas);
    etas = max(0,etas);
    lambdas = 1./etas;  % Just to make sure...
    phi = phi./sqrt(delta); % This removes the impact of delta (w).
    
%% For debugging
%     N2 = length(t);
%     phiX = zeros([N2,L]);
%     for l=1:L
%         for n=1:N2
%             phiX(n,l) = lambdas(l).*(delta.*K(n,:)*phi(:,l)); 
%         end
%     end
% 
%     disp('test');
end

