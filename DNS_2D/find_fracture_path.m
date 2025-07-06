function [bestPathDetails] = find_fracture_path(a,b,system_rotation_angle,X0)

crack_sp = a; F = b;
F = F - crack_sp;
crack_sp = [0 0] ;
x1 = a(1);  y1 = a(2);

centroid = mean(X0);
sorting_angles = atan2(X0(:,2) - centroid(2), X0(:,1) - centroid(1));
[~, order] = sort(sorting_angles);
X0 = X0(order, :);X0_closed = [X0; X0(1,:)];

if system_rotation_angle > 90
    system_rotation_angle = 90 - system_rotation_angle;
end

[K_IC,rotated_planes,planes,sortedIndices,SIF,path,allBasis] = basis_finder(system_rotation_angle);

paths = rotated_planes(sortedIndices,:);
allBasis = allBasis(sortedIndices,:);
rotated_planes = rotated_planes(sortedIndices,:);

allPaths = {};  planeContributions = [];   effectiveKICs = [];
applied_stress_direction = [0 100];
angle_in_between = acosd(dot(applied_stress_direction,F)/(norm(F)*(norm(applied_stress_direction))));

if angle_in_between == 0 || angle_in_between == 180
    %     disp('No Valid Path -x1');
    %     disp(angle_in_between);
    allPaths = {};
else

    for indexy = 1:length(sortedIndices)
        for j = 1:length(sortedIndices)
            combinedBasis = [allBasis{(indexy)}; allBasis{(j)}];

            if rank(combinedBasis) == 2
                coeffs = linsolve(combinedBasis', F');

                if any(coeffs < -10e-15)
                    continue;
                else
                    v1 = combinedBasis(1,:) * coeffs(1);
                    v2 = combinedBasis(2,:) * coeffs(2);
                    combinedKIC = SIF([indexy,j], :);
                    contributions = abs(coeffs) / sum(abs(coeffs)) * 100;
                    Plane_C1 = sum(contributions(1));
                    Plane_C2 = sum(contributions(2));
                    currentPathInfo = struct('planeIndices', [sortedIndices(indexy), sortedIndices(j)], ...
                        'planes', {rotated_planes([indexy,j], :)},'parent',{planes([sortedIndices(indexy), sortedIndices(j)], :)}, ...
                        'KIC', combinedKIC, 'basis', combinedBasis, 'coeffs', coeffs, 'contributions', contributions, 'planeC1', Plane_C1, 'planeC2', Plane_C2);
                    allPaths{end+1} = currentPathInfo;
                end
            end
        end
    end
end

pathsToRemove = [];

for pathIdx = 1:length(allPaths)
    currentPath = allPaths{pathIdx};
    v1 = currentPath.basis(1,:) * currentPath.coeffs(1);
    v2 = currentPath.basis(2,:) * currentPath.coeffs(2);

    endPoint1 = [v1(1)+x1, v1(2)+y1];
    endPoint2 = [v1(1)+v2(1)+x1, v1(2)+v2(2)+y1];
    inside1 = inpolygon(endPoint1(1), endPoint1(2), X0_closed(:,1), X0_closed(:,2));
    inside2 = inpolygon(endPoint2(1), endPoint2(2), X0_closed(:,1), X0_closed(:,2));

    if inside1 && inside2
        KIC1 = K_IC(currentPath.planeIndices(1));
        KIC2 = K_IC(currentPath.planeIndices(2));
        effectiveKIC = (currentPath.planeC1 / 100) * KIC1 + (currentPath.planeC2 / 100) * KIC2;
        allPaths{pathIdx}.Eff_KIC = effectiveKIC;
        effectiveKICs(end+1) = effectiveKIC;
    else
        pathsToRemove = [pathsToRemove, pathIdx];
    end
end

allPaths(pathsToRemove) = [];
[lowestEffectiveKIC, iidx] = min(effectiveKICs);

if ~isempty(allPaths)
    bestPath = allPaths{iidx};
    v1 = bestPath.basis(1,:) * bestPath.coeffs(1);
    v2 = bestPath.basis(2,:) * bestPath.coeffs(2);
else
    bestPathDetails = [];
    return;
end

bestPathDetails = struct();
bestPathDetails.S = a;
bestPathDetails.F = F+a;
bestPathDetails.I = [v1(1)+x1 v1(2)+y1];
bestPathDetails.SIF = K_IC(bestPath.planeIndices);
bestPathDetails.Effective_KIC = lowestEffectiveKIC;
bestPathDetails.planes = bestPath.planes;
bestPathDetails.parent = bestPath.parent;
bestPathDetails.contributions = [bestPath.planeC1, bestPath.planeC2];

end