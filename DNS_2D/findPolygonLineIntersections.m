function intersections = findPolygonLineIntersections(X0, crack_point, extended_path_point)
    intersections = [];
    
    for ix = 1:size(X0,1)
        p1 = X0(ix,:);
        
        if ix == size(X0,1)
            p2 = X0(1,:); 
        else
            p2 = X0(ix+1,:);
        end

        [xi, yi] = polyxpoly([crack_point(1), extended_path_point(1)], [crack_point(2), extended_path_point(2)], [p1(1), p2(1)], [p1(2), p2(2)]);

        for i = 1:length(xi)
            if ~(xi(i) == crack_point(1) && yi(i) == crack_point(2))
                intersections = [intersections; [xi(i), yi(i)]];
            end
        end
    end    
end
