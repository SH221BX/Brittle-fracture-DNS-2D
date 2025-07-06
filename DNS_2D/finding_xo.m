function [XX] = finding_xo(V, C, affected_cell, bx, constraints)

cell_verts = V(C{affected_cell}, :);
cell_verts = cell_verts(~isinf(cell_verts(:,1)), :);

corners = [0 0; bx 0; bx bx; 0 bx];
bufferZone = 1e-3;

for vv = 1:size(corners, 5)
    if shouldIncludeCorner(cell_verts, corners(vv,:), bufferZone)
        cell_verts = [cell_verts; corners(vv,:)];
    end
end

if inpolygon(0, 0, cell_verts(:,1), cell_verts(:,2))
    cell_verts = [cell_verts;[0,0]];
elseif inpolygon(bx, bx, cell_verts(:,1), cell_verts(:,2))
    cell_verts = [cell_verts;[bx,bx]];
elseif inpolygon(0, bx, cell_verts(:,1), cell_verts(:,2))
    cell_verts = [cell_verts;[0,bx]];
elseif inpolygon(bx, 0, cell_verts(:,1), cell_verts(:,2))
    cell_verts = [cell_verts;[bx,0]];
end

inBox_indices = [];
inBox_indices = (cell_verts(:,1) >= 0) & (cell_verts(:,1) <= bx) & (cell_verts(:,2) >= 0) & (cell_verts(:,2) <= bx);
XX = cell_verts(inBox_indices, :);

cell_int_points = [];
numVertices = size(cell_verts, 1);

for j = 1:numVertices - 1
    [xi, yi] = polyxpoly(cell_verts(j:j+1, 1), cell_verts(j:j+1, 2), constraints(:,1), constraints(:,2));
    if ~isempty(xi) && ~isempty(yi)
        cell_int_points = [cell_int_points; xi, yi];
    end
end

if size(cell_verts, 1) > 1
    [xi, yi] = polyxpoly([cell_verts(end, 1), cell_verts(1, 1)], [cell_verts(end, 2), cell_verts(1, 2)], constraints(:,1), constraints(:,2));
    if ~isempty(xi) && ~isempty(yi)
        cell_int_points = [cell_int_points; xi, yi];
    end
end


XX = [XX; cell_int_points];
XX = unique(XX,'rows','stable');

end

function isCornerIncluded = shouldIncludeCorner(cellVerts, corner, bufferZone)
cornerBuffer = [corner(1) - bufferZone, corner(2) - bufferZone;
    corner(1) + bufferZone, corner(2) - bufferZone;
    corner(1) + bufferZone, corner(2) + bufferZone;
    corner(1) - bufferZone, corner(2) + bufferZone;
    corner(1) - bufferZone, corner(2) - bufferZone];

isCornerIncluded = false;
for i = 1:size(cellVerts, 1)-1
    edgeStart = cellVerts(i, :);
    edgeEnd = cellVerts(i+1, :);

    for j = 1:4
        bufferStart = cornerBuffer(j, :);
        bufferEnd = cornerBuffer(j+1, :);

        [xi, yi] = polyxpoly([edgeStart(1), edgeEnd(1)], [edgeStart(2), edgeEnd(2)],...
            [bufferStart(1), bufferEnd(1)], [bufferStart(2), bufferEnd(2)]);
        if ~isempty(xi)
            isCornerIncluded = true;
            return;
        end
    end
end
end

