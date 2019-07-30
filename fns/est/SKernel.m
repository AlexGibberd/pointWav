function [Shat,Chat] = SKernel(E,A,b,wavType,kappa,SupK)
%SKERNEL Spectrum estimation via direct sampling

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
if(P>1)
    % Find events in approximate support and scale
    for i=1:P
        for j=1:P
            E_shift{i} = shiftScale(E{i},A,B,SupK);
            E_shift{j} = shiftScale(E{j},A,B,SupK);
        end
    end


    for i=1:P
        for j=1:P
            Ei = E_shift{i};
            Ej = E_shift{j};
            [EI,EJ] = meshgrid(Ei,Ej);
            if(strcmp(wavType,'Morlet')) 
                % We assume d=1..
                K = errFnKernelmorlet(EI,EJ,kappa);
            elseif(strcmp(wavType,'Mexican'))
                % We assume scale = 1, sigma = 1
                K = mexicanKernel(EI,EJ,1,1,kappa);
            end

            % Scale sum according to A
            % kernels are defined for A=1
            Shat(i,j) = sum(sum(K))/A;
        end
    end
    for i=1:P
       for j=1:P
           Chat(i,j) = abs(Shat(i,j))^2./(Shat(i,i)*Shat(j,j)) ;
       end
    end
    
else
    E_shift = shiftScale(E,A,B,SupK);
    [EI,EJ] = meshgrid(E_shift,E_shift);
    if(strcmp(wavType,'Morlet')) 
        % We assume d=1..
        K = errFnKernelmorlet(EI,EJ,kappa);
    elseif(strcmp(wavType,'Mexican'))
        % We assume scale = 1, sigma = 1
        K = mexicanKernel(EI,EJ,1,1,kappa);
    end

    % Scale sum according to A
    % kernels are defined for A=1
    Shat = sum(sum(K))/A;
    Chat = 0;
end


end

function [xp] = shiftScale(x,A,B,SupK)
% Evaluate the grid and sample kernel
    xp = (x-B)./A; % shift and scale data
    % Filter data to approximate support if given
    if(SupK>0)
        filter = logical(abs(xp)>SupK);
        x_start = find(filter==0,1,'first');
        x_end = find(filter==0,1,'last');
        xp(filter) = [];
        Nx = length(xp);
    end
end