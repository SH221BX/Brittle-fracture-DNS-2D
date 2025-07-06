function [affectedcell1, affectedcell2, points_in_between] = findSharedGrainsByEdge(givenPoint, V, C, cell_orientations, bx)
    epsilon = 1e-4;  % Precision for geometric calculations
    constraints = [0 0; 0 1; 1 1; 1 0; 0 0] * bx;
    affectedcell1 = -1;
    affectedcell2 = -1;
    points_in_between = [];
    sharedGrains = [];

    for i = 1:numel(C)
        vertices = V(C{i}, :);
        vertices = vertices(~isinf(vertices(:,1)), :);
        for j = 1:size(vertices, 1)
            nextj = mod(j, size(vertices, 1)) + 1;
            edgeStart = vertices(j, :);
            edgeEnd = vertices(nextj, :);

            if isPointOnLine(givenPoint, edgeStart, edgeEnd, epsilon)
                sharedGrains = unique([sharedGrains; i]);  
                % Store unique grain indices
                points_in_between = [points_in_between; edgeStart; edgeEnd];  
                % Collect edge points
            end
        end
    end

    if numel(sharedGrains) == 1
        affectedcell1 = sharedGrains(1);
        affectedcell2 = -1;  
        % No second cell affected
    elseif numel(sharedGrains) >= 2
        affectedcell1 = sharedGrains(1);
        affectedcell2 = sharedGrains(2);
    end
end


function onLine = isPointOnLine(point, lineStart, lineEnd, epsilon)
    v = lineEnd - lineStart;
    w = point - lineStart;

    c1 = dot(w,v);
    if c1 <= 0
        distance = norm(point - lineStart);
    else
        c2 = dot(v,v);
        if c2 <= c1
            distance = norm(point - lineEnd);
        else
            b = c1 / c2;
            Pb = lineStart + b * v;
            distance = norm(point - Pb);
        end
    end

    onLine = distance < epsilon;
end
