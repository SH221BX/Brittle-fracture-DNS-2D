function shared_cells = find_Shared_grains(point, V, C)
    tolerance = 1e-3;  
    % Adjusted tolerance for geometric calculations
    shared_cells = [];

    % Iterate through each cell
    for p = 1:length(C)
        cell_vertices = V(C{p}, :);
        n_vertices = size(cell_vertices, 1);

        % Check each edge of the cell
        for j = 1:n_vertices
            % Indices for current and next vertex to form an edge
            current_vertex = cell_vertices(j, :);
            next_vertex = cell_vertices(mod(j, n_vertices) + 1, :);

            % Check if the point is on the line segment
            if isPointOnLineSegment(point, current_vertex, next_vertex, tolerance)
                shared_cells = [shared_cells, p];
                break;  % Stop checking further edges if one match is found
            end
        end
    end
end

function onLine = isPointOnLineSegment(point, v1, v2, tolerance)
    lineVec = v2 - v1;
    pointVec = point - v1;

    % Line vector length squared (avoiding sqrt for performance)
    lineLenSquared = dot(lineVec, lineVec);

    % Projection of pointVec onto lineVec
    projection = dot(pointVec, lineVec);

    % Check if the projection lies within the line segment
    if projection < 0 || projection > lineLenSquared
        onLine = false;
        return;
    end

    % Distance from the point to the line (using cross product magnitude)
    distance = norm(cross([lineVec, 0], [pointVec, 0])) / sqrt(lineLenSquared);

    % Check distance against tolerance
    onLine = (distance <= tolerance);
end
