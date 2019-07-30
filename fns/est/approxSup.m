function [S] = approxSup(wavType,alpha,kappa)
    % S - approx half width of kernel (assumes symmetrical about s=t)
    % Assumes kernel K() centered on zero
    N = 10000; % Big number of samples
    lim = 100; % Maximum range of support to probe
    Tl = linspace(0,lim,N);

    if(strcmp(wavType,'Morlet')) 
        % We assume d=1..
        for i=1:N
            Kd(i) = errFnKernelmorlet(Tl(i),Tl(i),kappa);
        end
    elseif(strcmp(wavType,'Mexican'))
        % We assume scale = 1, sigma = 1
        for i=1:N
            Kd(i) = mexicanKernel(Tl(i),Tl(i),1,1,kappa);
        end
        % For large t, this can be nan
        Kd(isnan(Kd))=0;
    else
        S = 0;
        return; 
    end
    
    K2 = Kd.*conj(Kd);
    % In below we assume wavelet is symmetric
    S = Tl(find( (cumsum(K2)./sum(K2)) > 1-(alpha/2), 1 ));
end
