clear all; format short;  figure; hold on;
tic;warning off;

for sd = 5:5
%     rng(15,'v4'); 
    GBC = 1;
    bx = 100; edges = [0 0 1 1 0;0 1 1 0 0] * bx;
    num_nodes = 50;
    num_nodes = num_nodes-4;
    nodes_pos = rand(num_nodes, 2)*bx;

    constraints=[0 0; 0 1; 1 1; 1 0;0 0] * bx;
    nodes_pos = [nodes_pos; constraints]; nodes_pos = unique(nodes_pos, 'rows');
    x = nodes_pos(:, 1); y = nodes_pos(:, 2);
    dt = delaunayTriangulation(x,y); [V, C] = voronoiDiagram(dt);

    for iteration = 1:40
        newX = zeros(size(x)); newY = zeros(size(y));
        for i = 1:length(C)
            if ~isempty(C{i})
                Vx = V(C{i}, 1);
                Vy = V(C{i}, 2);
                Xa = [Vx(2:end);Vx(1)]; Ya = [Vy(2:end);Vy(1)];
                A = 1/2*sum(Vx.*Ya-Xa.*Vy);
                cx = (1/(6*A)*sum((Vx + Xa).*(Vx.*Ya-Xa.*Vy)));
                cy = (1/(6*A)*sum((Vy + Ya).*(Vx.*Ya-Xa.*Vy)));
                if inpolygon(cx, cy, constraints(:,1), constraints(:,2))
                    newX(i) = cx;
                    newY(i) = cy;
                else
                    newX(i) = x(i);
                    newY(i) = y(i);
                end
            end
        end
        x = newX; y = newY;
        [V, C] = VoronoiBounded(x, y, constraints);
    end
  
    cell_orientations = -90 + (180) * rand(1, size(C, 1));
    all_cell_verts = []; int_points = [];grains_lib = []; colormap('hsv');

    for cell_no = 1:size(C, 1)
        [cell_verts] = finding_xo(V, C, cell_no, bx, constraints);
        centroid = mean(cell_verts);
        sorting_angles1 = atan2(cell_verts(:,2) - centroid(2), cell_verts(:,1) - centroid(1));
        [~, order1] = sort(sorting_angles1); cell_verts = cell_verts(order1, :);
        all_cell_verts = [all_cell_verts;cell_verts];
          orientation = cell_orientations(cell_no);
          color_value = orientation;
        patch(cell_verts(:,1), cell_verts(:,2),color_value,'FaceAlpha',0, 'EdgeColor','k');hold on;
        grain_areas(cell_no) = polyarea(cell_verts(:,1), cell_verts(:,2));
        grain_diameter(cell_no)=sqrt((4*grain_areas(cell_no))/(3.1416));
        grains_lib(cell_no).index = cell_no;
        grains_lib(cell_no).vertices = cell_verts;
        grains_lib(cell_no).diameter = grain_diameter(cell_no);
        grains_lib(cell_no).area = grain_areas(cell_no);

        plot(cell_verts(:,1), cell_verts(:,2), '-k');hold on;

                Cx = mean(cell_verts(:,1)); Cy = mean(cell_verts(:,2));
                if (Cx >= 0 && Cx <= bx) && (Cy >= 0 && Cy <= bx)
                    orientationText = sprintf('%.2f°', orientation);
                 end
    end


    grain_diameter = [grains_lib.diameter]';
    grain_area = [grains_lib.area]';

    Average_area = 10000/length(grain_diameter);
    Average_D = sqrt(4*Average_area/3.1416);  % disp(Average_D);


    crack_point = [0, bx/2];
    scatter(crack_point(1),crack_point(2),30,'s', 'filled', 'MarkerEdgeColor', 'k','MarkerFaceColor','red');hold on;

    colors = [[1 1 0]; [1 0 0]; [0 0 1]; [0.5 0 0.5]];
    all_paths = [];  n=0; GB_path = []; Grain_path = []; all_GB_paths = [];
    local_GB_paths = []; local_GB_path = []; temp_path_option = [];
    local_T_path = []; local_TGB_path = []; TGB_all = []; TGB_path = [];
    all_TGB_paths = []; T_all_path = []; cracked_grain = []; crack_path = [];

     while true
                all_paths = [];   GB_path = []; Grain_path = []; all_GB_paths = [];
                local_GB_paths = []; local_GB_path = []; temp_path_option = [];
                local_T_path = []; local_TGB_path = []; TGB_all = []; TGB_path = [];
                all_TGB_paths = []; T_all_path = [];
                V_rounded = round(V, 6); crack_point_rounded = round(crack_point, 6);
                m = ismember(crack_point_rounded, V_rounded, 'rows');

                if any(m)
                    shared_grains = find_Shared_grains(crack_point, V, C);

                    % To next triple junction
                    % _______________________

                    for grain_number=1:length(shared_grains)
                        [X0] = finding_xo(V, C, shared_grains(grain_number), bx, constraints);
                        centroid = mean(X0);
                        sorting_angles = atan2(X0(:,2) - centroid(2), X0(:,1) - centroid(1));
                        [~, order] = sort(sorting_angles);
                        X0 = X0(order, :);X0_closed = [X0; X0(1,:)];
                        index_of_crack_point = find(ismember(round(X0, 6), crack_point_rounded, 'rows'));

                        if index_of_crack_point == 1
                            connected_points = X0([end, 2], :);
                        elseif index_of_crack_point == size(X0, 1)
                            connected_points = X0([index_of_crack_point-1, 1], :);
                        else
                            connected_points = X0([index_of_crack_point-1, index_of_crack_point+1], :);
                        end

                        connected_points_rounded = round(connected_points, 8);
                        cdx = find(connected_points_rounded(:,1) <= crack_point(1), 2);
                        connected_points(cdx,:)=[];

                        if isempty(connected_points)
                            % disp('No connected points');
                            continue;

                        elseif ~isempty(connected_points)
                            cp = 0;
                            for cp = 1:size(connected_points,1)
                                % Destination: triple junction - following GB
                                % -----------------------------------------
                                local_GB_path.S = crack_point;
                                local_GB_path.F = connected_points(cp,:);
                                local_GB_path.I = (local_GB_path.S+local_GB_path.F)/2;
                                [affectedcell1, affectedcell2,points_in_between] = findSharedGrainsByEdge(crack_point, V, C,cell_orientations,bx);
                                local_GB_vector = connected_points(cp,:) - crack_point;
                                local_angle_wrt_x = rad2deg(acos(local_GB_vector(1)/norm(local_GB_vector)));
                                local_GB_path.SIF = GBC/(((cosd(local_angle_wrt_x)).^1));
                                local_GB_path.Effective_KIC = GBC/(((cosd(local_angle_wrt_x)).^1));
                                local_GB_path.planes = local_GB_vector;
                                local_GB_path.parent = local_GB_vector;
                                local_GB_path.contributions = 100;
                                local_GB_path.grains = [affectedcell1, affectedcell2];   %shared_grains(grain_number);
                                temp_path_option = [temp_path_option;local_GB_path];
                     
                                % Destination: Next triple junction - using planes in the grain
                                % -------------------------------------------------------------
                                theta = cell_orientations(shared_grains(grain_number));
                                transgranular_path = find_fracture_path(crack_point,connected_points(cp,:),theta,X0);
                                if ~isempty(transgranular_path)
                                    local_TGB_path = transgranular_path;
                                    local_TGB_path.grains = [shared_grains(grain_number)];
                                    temp_path_option = [temp_path_option;local_TGB_path];
                                else
                                    % fprintf('No T path to NTJ\n');
                                    continue;
                                end
                            end
                        end

                        % Source: triple junction - Transgranular paths
                        % ---------------------------------------------

                        theta = cell_orientations(shared_grains(grain_number));
                        [K_IC,rotated_planes,planes,sortedIndices,SIF,path,allBasis] = basis_finder(theta);
                        paths = bx*rotated_planes(sortedIndices,:);

                        for loop = 1:length(paths)
                            extended_path_point = [crack_point(1)+paths(loop,1)   crack_point(2)+paths(loop,2)];
                            if extended_path_point(1)>=crack_point(1)
                                intersections = findPolygonLineIntersections(X0, crack_point, extended_path_point);
                                for i = size(intersections, 1):-1:1
                                    if abs(intersections(i, 1) - crack_point(1)) < (1e-4) && abs(intersections(i, 2) - crack_point(2)) < (1e-4)
                                        intersections(i, :) = [];
                                    end
                                end

                                if isempty(intersections)
                                    continue;

                                elseif (~isempty(intersections)) && (intersections(1)>=crack_point(1))
                                    local_T_path.S = crack_point;
                                    local_T_path.F = intersections;
                                    local_T_path.I = (local_T_path.S+local_T_path.F)/2;
                                    local_T_vector = intersections - crack_point;
                                    local_T_path.SIF = SIF(loop);
                                    local_T_path.Effective_KIC = SIF(loop);
                                    local_T_path.planes = rotated_planes(loop,:);
                                    local_T_path.parent = planes([sortedIndices(loop)],:);
                                    local_T_path.contributions = 100;
                                    local_T_path.grains = shared_grains(grain_number);
                                    temp_path_option = [temp_path_option;local_T_path];
                                end
                            end
                        end
                    end

                    % temp_path_option = makeStructUnique(temp_path_option);
                    i=0; toRemoveIdx = [];

                    for i = 1:length(temp_path_option )
                        if temp_path_option(i).F(1) <= crack_point(1)
                            toRemoveIdx = [toRemoveIdx, i];
                        end
                    end
                    if ~isempty(toRemoveIdx)
                        temp_path_option(toRemoveIdx) = [];
                    end

                    E_KIC_values = arrayfun(@(x) x.Effective_KIC, temp_path_option);
                    [~, SIndices] = sort(E_KIC_values);
                    temp_path_option = temp_path_option(SIndices);

                    The_path = temp_path_option(1);

                    if length(The_path.grains)<2
                        cracked_grain = [cracked_grain,The_path.grains];
                        % plot([The_path.S(1),The_path.F(1)],[The_path.S(2),The_path.F(2)],':','Color','r','LineWidth',2.5); hold on
                        % scatter(The_path.F(1), The_path.F(2), 40, 's', 'filled', 'MarkerEdgeColor', 'k','MarkerFaceColor','blue');hold on;
                    else
                        % plot([The_path.S(1),The_path.F(1)],[The_path.S(2),The_path.F(2)],':','Color','g','LineWidth',2.5); hold on
                        % scatter(The_path.F(1), The_path.F(2), 40, 's', 'filled', 'MarkerEdgeColor', 'k','MarkerFaceColor','blue');hold on;
                    end
%                     scatter(crack_point(1),crack_point(2),50,'filled','ko');hold on;
                    crack_point = The_path(1).F;
                    %  fprintf('crack_point %.4f %.4f \n',crack_point(1),crack_point(2));
                    

                else

                    % Not a triple junction
                    % ---------------------

                    affectedcell2 = -1; affectedcell1 =-1;
                    % disp('Not a triple junction');
                    [affectedcell1, affectedcell2,points_in_between] = findSharedGrainsByEdge(crack_point, V, C,cell_orientations,bx);
                    affected_cells = [affectedcell1, affectedcell2];
                    mdx = find(affected_cells(:)>0);
                    idx = find(round(points_in_between(:,1),6) >= round(crack_point(1),6), 2);

                    % Grain boundary fracture paths
                    % ------------------------------

                    for z1 = 1:length(idx)

                        point_of_interest = points_in_between(idx(z1), :);

                        xint = 0; yint = 0;
                        if any(point_of_interest < 0) || any(point_of_interest > bx)
                            [xint, yint] = polyxpoly([crack_point(1), point_of_interest(1)], [crack_point(2), point_of_interest(2)], ...
                                constraints(:,1), constraints(:,2));
                            if ~isempty(xint) && ~isempty(yint)
                                point_of_interest = [xint(1), yint(1)];
                            else
                                disp('Error: No intersection found or outside boundary.');
                            end
                        end

                        GB_all.S = crack_point;
                        GB_all.F = point_of_interest;
                        GB_all.I = (GB_all.S+GB_all.F)/2;
                        GB_vector = GB_all.F - GB_all.S;
                        local_angle_wrt_x = rad2deg(acos(GB_vector(1)/norm(GB_vector)));
                        GB_all.SIF = GBC/(((cosd(local_angle_wrt_x)).^1));
                        GB_all.Effective_KIC = GBC/(((cosd(local_angle_wrt_x)).^1));
                        GB_all.planes = GB_vector;
                        GB_all.parent = GB_vector;
                        GB_all.contributions = 100;
                        GB_all.grains = affected_cells(affected_cells(:)>0);
                        TGB_all = [TGB_all;GB_all];

                           
                        for Q = 1:length(mdx)
                            [X0] = finding_xo(V, C, affected_cells(mdx(Q)), bx, constraints);
                            centroid = mean(X0); sorting_angles = atan2(X0(:,2) - centroid(2), X0(:,1) - centroid(1));
                            [~, order] = sort(sorting_angles); X0 = X0(order, :);X0_closed = [X0; X0(1,:)];
                            theta = cell_orientations(affected_cells(mdx(Q)));

                            TGB_path_TX= find_fracture_path(crack_point,points_in_between(idx(z1), :),theta,X0);
                            if ~isempty(TGB_path_TX)
                                TGB_path_TX.grains = [affected_cells(mdx(Q))];
                                TGB_all = [TGB_all;TGB_path_TX];
                            end
                        end
                    end


                    % Transgranular fracture paths
                    % ----------------------------
                    for q = 1:length(mdx)

                        if ~ismember(affected_cells(mdx(q)), cracked_grain)
                            [X0] = finding_xo(V, C, affected_cells(mdx(q)), bx, constraints);

                            centroid = mean(X0);
                            sorting_angles = atan2(X0(:,2) - centroid(2), X0(:,1) - centroid(1));
                            [~, order] = sort(sorting_angles);
                            X0 = X0(order, :);X0_closed = [X0; X0(1,:)];

                            theta = cell_orientations(affected_cells(mdx(q)));
                            [K_IC,rotated_planes,planes,sortedIndices,SIF,path,allBasis] = basis_finder(theta);
                            paths = bx*rotated_planes(sortedIndices,:);


                            for loop = 1:length(paths)
                                extended_path_point = [crack_point(1)+paths(loop,1)   crack_point(2)+paths(loop,2)];
                                intersections = findPolygonLineIntersections(X0, crack_point, extended_path_point);

                                for i = size(intersections, 1):-1:1
                                    if abs(intersections(i, 1) - crack_point(1)) < (1e-3) && abs(intersections(i, 2) - crack_point(2)) < (1e-3)
                                        intersections(i, :) = [];
                                    end
                                end

                                if (~isempty(intersections)) && (intersections(1)>=crack_point(1))
                                    T_path.S = crack_point;
                                    T_path.F = intersections;
                                    T_path.I = (T_path.S+T_path.F)/2;
                                    T_vector = intersections - crack_point;
                                    T_path.SIF = SIF(loop);
                                    T_path.Effective_KIC = T_path.SIF;
                                    T_path.planes = rotated_planes(loop,:);
                                    T_path.parent = planes([sortedIndices(loop)],:);
                                    T_path.contributions = 100;
                                    T_path.grains = affected_cells(mdx(q));
                                    T_all_path = [T_all_path;T_path];

                                    toRemoveIdx = [];i=0;
                                    for i = 1:length(T_all_path)
                                        if T_all_path(i).F(1) <= crack_point(1)
                                            toRemoveIdx = [toRemoveIdx, i];
                                        end
                                    end
                                    if ~isempty(toRemoveIdx)
                                        T_all_path(toRemoveIdx) = [];
                                    end        

                                end

                                if isempty(T_all_path)
                                    continue;
                                end
                            end

                            TGB_all = [TGB_all;T_all_path];

                            TGB_all= makeStructUnique(TGB_all);
                            toRemoveIdx = [];
                            i=0;
                            for i = 1:length(TGB_all)
                                if round(TGB_all(i).F(1),6) <= round(crack_point(1),6)
                                    toRemoveIdx = [toRemoveIdx, i];
                                end
                            end

                            if ~isempty(toRemoveIdx)
                                TGB_all(toRemoveIdx) = [];
                            end

                            E_SIF_values = arrayfun(@(x) x.Effective_KIC, TGB_all);
                            [~, SI] = sort(E_SIF_values);
                            TGB_all= TGB_all(SI);
                        end
                    end

                    The_path = TGB_all(1);
                    cracked_grain = [cracked_grain,The_path.grains];
                    crack_point = The_path(1).F;

                end
                if length(The_path.grains)>1
                     plot([The_path.S(1),The_path.F(1)],[The_path.S(2),The_path.F(2)],'-','Color','g','LineWidth',2,'LineStyle','-','AlignVertexCenters','on'); hold on
                else
                     plot([The_path.S(1),The_path.I(1)],[The_path.S(2),The_path.I(2)],'-','Color','r','LineWidth',2); hold on
                     plot([The_path.I(1),The_path.F(1)],[The_path.I(2),The_path.F(2)],'-','Color','r','LineWidth',2); hold on
                end
                crack_path = [crack_path;The_path];

                if round(crack_point(1),4) >= bx
                    break;
                end
     end
end
            scatter(The_path.F(1), The_path.F(2), 30, 's', 'filled', 'MarkerEdgeColor', 'k','MarkerFaceColor','blue');hold on;

            plot(edges(1,:), edges(2,:), 'k','LineWidth', 0.1); hold on;
            xticks([]);yticks([]);
            xlabel('');ylabel('');
            axis equal;
            grid off;box on; xlim([0 bx]); ylim([0 bx]);
            set(gcf, 'Renderer', 'opengl');

elapsedTime = toc;
disp(['Elapsed time: ', num2str(elapsedTime), ' seconds.']);


