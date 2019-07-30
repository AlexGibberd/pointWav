function [phiX] = eigenWav(phi,lambdas,s,x,kappa,a,b,varargin)
%EIGENWAV evaluates samples from the Nystrom extension
%   
% Required: 
%   psi (Ns xL) vector of Nystrom sampled eigen-wavelets (or single value containing sum)
%   lambdas - L vector with eigenvalues
%   s - nystrom sample points
%   x - evaluation points
%   kappa - smoothing size
% Optional:
%   w - weight vector (of same length as s), default, derives this from s
%   wavType - 'Morlet' (default), or 'Mexican'
%   SupK - Approximate (half) support of the kernel
%   padding - 1/0 (default 1) only return samples in support, or pad to
%   original size
%   transform - 1/0 (default 0), whether to return sum, i.e. eigenWav transform,
%   or evaluation
%
% Output:
%   phiX (Nx x L)

% First parse inputs. Make default weight vector
    Ns = length(s);
    Nx = length(x);
    L = length(lambdas);
    default_w = ((max(s)-min(s))/Ns).*ones([1,Ns]);
    pa = inputParser;  
    addOptional(pa,'w',default_w,@isnumeric);
    addOptional(pa,'wavType','Morlet');
    addOptional(pa,'SupK',0,@isnumeric);
    addOptional(pa, 'padding', 1, @isnumeric);
    addOptional(pa, 'transform', 0, @isnumeric);

    parse(pa,varargin{:});  
    w = pa.Results.w;
    wavType = pa.Results.wavType;
    SupK = pa.Results.SupK;
    padding = pa.Results.padding;
    transform = pa.Results.transform;
    
% Evaluate the grid and sample kernel
    xp = (x-b)./a; % shift and scale data
% Filter data to approximate support if given
    if(SupK>0)
        filter = logical(abs(xp)>SupK);
        x_start = find(filter==0,1,'first');
        x_end = find(filter==0,1,'last');
        xp(filter) = [];
        Nx = length(xp);
    end
    [X,S] = meshgrid(xp,s);
    if(strcmp(wavType,'Morlet')) 
        % We assume d=1..
        Kx = errFnKernelmorlet(X,S,kappa);        
    elseif(strcmp(wavType,'Mexican'))
        % We assume scale = 1, sigma = 1
        Kx = mexicanKernel(X,S,1,1,kappa);
    end
% Perform weighting of kernel
    if(transform == 0)
        % Evaluate eigenwav at points
        phiX = zeros([Nx,L]);
        for l=1:L
            for n=1:Nx
                phiX(n,l) = lambdas(l)*(diag(w)*Kx(:,n))'*phi(:,l); 
            end
        end
        phiX = phiX./sqrt(a); % Rescale output according to scale relation
    else
        % Take transform using eigenwav
        for l=1:L
            phiX(l) = sum(lambdas(l)*(diag(w)*Kx)'*phi(:,l));
        end
        phiX = phiX./sqrt(a);
        return;
    end

% If we trimmed the data before, we need to add zeros back
    if(SupK>0 && padding == 1)
        phiX2 = zeros(length(x),L);
        phiX2(x_start:x_end,:) = phiX;
        phiX = phiX2;
    end
end

