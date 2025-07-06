function [intersect, point] = lineIntersect(P1, P2, Q1, Q2)
    % Initialize
    intersect = false;
    point = [NaN, NaN];
    
    tol = 1e-9;  % Define a tolerance

    % Line P represented as P1 + t(P2 - P1)
    % Line Q represented as Q1 + u(Q2 - Q1)
    A1 = P2(2) - P1(2);
    B1 = P1(1) - P2(1);
    C1 = A1 * P1(1) + B1 * P1(2);

    A2 = Q2(2) - Q1(2);
    B2 = Q1(1) - Q2(1);
    C2 = A2 * Q1(1) + B2 * Q1(2);

    det = A1 * B2 - A2 * B1;

    if abs(det) < tol
        return; % lines are parallel (or coincident)
    end
    
    x = (B2 * C1 - B1 * C2) / det;
    y = (A1 * C2 - A2 * C1) / det;

    % Check if the point is on both line segments with a tolerance
    if (x >= min(P1(1), P2(1)) - tol && x <= max(P1(1), P2(1)) + tol) && ...
       (y >= min(P1(2), P2(2)) - tol && y <= max(P1(2), P2(2)) + tol) && ...
       (x >= min(Q1(1), Q2(1)) - tol && x <= max(Q1(1), Q2(1)) + tol) && ...
       (y >= min(Q1(2), Q2(2)) - tol && y <= max(Q1(2), Q2(2)) + tol)
        intersect = true;
        point = [x, y];
    end
end

