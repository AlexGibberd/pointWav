function [Z] = constructZ(a,Amax,supK,T,oversamp)
    Na = length(a);
    for j=1:Na
        A(j) = a(j)*Amax;
        Sup(j) = psiSup(a(j),Amax,supK); % 1/2 Support at a(j) scale
    end

    Z = cell([Na,1]);   % Define scale/time sampling position
    for j=1:Na
        % Calculate number of samples to take at each scale a
        Nz = ceil((T-2*Sup(j))/(2*Sup(j)));     % Number of samples at scale a(j)
        % We allow sampling points to overlap inclusive of smoothing but
        % not in terms of wavelet support...
        Nz = max(oversamp*Nz,1); % Oversample...
        Z{j} = linspace(Sup(j)/T,1-(Sup(j)/T),Nz); % Assume symmetric with b=0 for now
        % Sample both sides of triangle
        if (Nz==1)
           Z{j}(2) = Z{j}(1);
           Z{j}(1) = 1-Z{j}(2);
        end
    end
    
end