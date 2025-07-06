function [X0, highlighted_vertices] = computeX0(V, C, affected_cell, bx, constraints, cell_orientations, current_crack_point)

    vertices = V(C{affected_cell}, :);
    vertices = vertices(~isinf(vertices(:,1)), :);  
    % disp(vertices);
    inBox_indices = [];
    inBox_indices = (vertices(:,1) >= 0) & (vertices(:,1) <= bx) & (vertices(:,2) >= 0) & (vertices(:,2) <= bx);
    XX = vertices(inBox_indices, :);
    cell_int_points = [];
    numVertices = size(vertices, 1);

    for j = 1:numVertices - 1
        [xi, yi] = polyxpoly(vertices(j:j+1, 1), vertices(j:j+1, 2), constraints(:,1), constraints(:,2));
        if ~isempty(xi) && ~isempty(yi)
            cell_int_points = [cell_int_points; xi, yi];
        end
    end

    if size(vertices, 1) > 1
        [xi, yi] = polyxpoly([vertices(end, 1), vertices(1, 1)], [vertices(end, 2), vertices(1, 2)], constraints(:,1), constraints(:,2));
        if ~isempty(xi) && ~isempty(yi)
            cell_int_points = [cell_int_points; xi, yi];
        end
    end

    highlighted_vertices = [];
    for i = 1:numVertices-1
        v1 = vertices(i, :);
        v2 = vertices(i + 1, :);
        if inpolygon(current_crack_point(1), current_crack_point(2), [v1(1), v2(1)], [v1(2), v2(2)])
            highlighted_vertices = [highlighted_vertices; v1; v2];
            break;
        end
    end
    
    if isempty(highlighted_vertices)
        v1 = vertices(numVertices, :);
        v2 = vertices(1, :);
        if inpolygon(current_crack_point(1), current_crack_point(2), [v1(1), v2(1)], [v1(2), v2(2)])
            highlighted_vertices = [highlighted_vertices; v1; v2];
        end
    end
    
    XX = [XX; cell_int_points];
    anglesd = rad2deg(atan2(XX(:,2) - current_crack_point(2), XX(:,1) - current_crack_point(1)));
    [~, sortedIndices] = sort(anglesd, 'descend','MissingPlacement','auto');
    X0 = XX(sortedIndices, :);
end
