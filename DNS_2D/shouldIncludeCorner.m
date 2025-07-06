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




