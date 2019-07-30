function [aidx, zidx] = findClosestSample(ESpec,coord)
    % Coord = [a,z]
    Na = length(ESpec.a);
    % Brute force search :/
    min_daz = 1000; % very lage..
    for j=1:Na
       for n=1:length(ESpec.Z{j})
           % Use l2 norm
           daz = sqrt((coord(1)-ESpec.a(j)).^2 + (coord(2)-ESpec.Z{j}(n)).^2);
           if daz<min_daz
              min_daz = daz;
              aidx = j;
              zidx = n;
           end
       end
    end
end
