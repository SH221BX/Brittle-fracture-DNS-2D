function [verticesOut, cellsOut] = VoronoiBounded(xPoints, yPoints, boundaryPoly)
% Compute bounding box
bounds = [min(xPoints), max(xPoints), min(yPoints), max(yPoints)]; %#ok<NASGU>

% Compute ranges and midpoint of the boundary polygon
rangeX = max(boundaryPoly(:,1)) - min(boundaryPoly(:,1));
rangeY = max(boundaryPoly(:,2)) - min(boundaryPoly(:,2));
rangeMax = max(rangeX, rangeY);
centerX = (max(boundaryPoly(:,1)) + min(boundaryPoly(:,1))) / 2;
centerY = (max(boundaryPoly(:,2)) + min(boundaryPoly(:,2))) / 2;

% Augment points to ensure bounded regions
augX = [xPoints; centerX + [0; 0; -5*rangeMax;  5*rangeMax]];
augY = [yPoints; centerY + [-5*rangeMax; 5*rangeMax; 0; 0]];

% Compute full Voronoi diagram
[vertices, cells] = voronoin([augX, augY]);

% Exclude the artificial points added
cellsOut = cells(1:end-4);
verticesOut = vertices;

% Clip each Voronoi cell by the boundary polygon
for cellIdx = 1:length(cellsOut)
    % Extract cell vertices and order clockwise
    [polyX, polyY] = poly2cw(verticesOut(cellsOut{cellIdx},1), verticesOut(cellsOut{cellIdx},2));
    % Intersection with boundary polygon
    [intersectX, intersectY] = polybool('intersection', boundaryPoly(:,1), boundaryPoly(:,2), polyX, polyY);

    % Find or append intersected vertices
    vertexIndices = nan(1, length(intersectX));
    for ptIdx = 1:length(intersectX)
        idxX = find(verticesOut(:,1) == intersectX(ptIdx));
        idxY = find(verticesOut(:,2) == intersectY(ptIdx));
        commonIdx = intersect(idxX, idxY);
        if ~isempty(commonIdx)
            vertexIndices(ptIdx) = commonIdx(1);
        else
            newIdx = size(verticesOut,1) + 1;
            verticesOut(newIdx,1) = intersectX(ptIdx);
            verticesOut(newIdx,2) = intersectY(ptIdx);
            vertexIndices(ptIdx) = newIdx;
        end
    end
    cellsOut{cellIdx} = vertexIndices;
end
end


