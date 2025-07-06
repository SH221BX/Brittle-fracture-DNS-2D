
clear all; clc;
aaa = []; bbb = []; ccc = [];

for seed_indice = 20:1:20

    rng(seed_indice,'v4'); GBC = 10000;
    bx = 100; edges = [0 0 1 1 0;0 1 1 0 0] *bx;
    num_nodes = 10000;
    nodes_pos = rand(num_nodes, 2)*bx;
    % num_nodes = num_nodes-4;
    % nodes_pos = unifrnd(0,1,2,num_nodes)*bx;
    % nodes_pos = nodes_pos';

    constraints=[0 0; 0 1; 1 1; 1 0;0 0] * bx;
    % nodes_pos = [nodes_pos; constraints];
    nodes_pos = unique(nodes_pos, 'rows');
    x = nodes_pos(:, 1); y = nodes_pos(:, 2);

    dt = delaunayTriangulation(x,y);
    %[V, C] = voronoiDiagram(dt);

    [V, C] = VoronoiBounded(x, y, constraints);

    % % figure; hold on;
    % convergenceThreshold = 0.1;
    % maxIterations = 100;
    % previousX = x; previousY = y;
    %
    % for iteration = 1:5
    %     newX = zeros(size(x));
    %     newY = zeros(size(y));
    %     for i = 1:length(C)
    %         if ~isempty(C{i})
    %             Vx = V(C{i}, 1);
    %             Vy = V(C{i}, 2);
    %             Xa = [Vx(2:end); Vx(1)];
    %             Ya = [Vy(2:end); Vy(1)];
    %             A = 1/2*sum(Vx.*Ya - Xa.*Vy);
    %             cx = (1/(6*A)*sum((Vx + Xa).*(Vx.*Ya - Xa.*Vy)));
    %             cy = (1/(6*A)*sum((Vy + Ya).*(Vx.*Ya - Xa.*Vy)));
    %             if inpolygon(cx, cy, constraints(:,1), constraints(:,2))
    %                 newX(i) = cx;
    %                 newY(i) = cy;
    %             else
    %                 newX(i) = x(i);
    %                 newY(i) = y(i);
    %             end
    %         end
    %     end
    %     maxMovement = max(sqrt((newX - previousX).^2 + (newY - previousY).^2));
    %     x = newX; y = newY; previousX = x;previousY = y;
    %     [V, C] = VoronoiBounded(x, y, constraints);
    % end


    figure; hold on; axis equal; box on;
    % plot(points(:,1), points(:,2), 'ro','MarkerSize',5);hold on;
    all_cell_verts = []; int_points = [];grains_lib = []; colormap('hsv');

    for cell_no = 1:size(C, 1)
        [cell_verts] = finding_xo(V, C, cell_no, bx, constraints);
        centroid = mean(cell_verts);
        sorting_angles1 = atan2(cell_verts(:,2) - centroid(2), cell_verts(:,1) - centroid(1));
        [~, order1] = sort(sorting_angles1); cell_verts = cell_verts(order1, :);
        all_cell_verts = [all_cell_verts;cell_verts];
        patch(cell_verts(:,1), cell_verts(:,2),'r','FaceAlpha',0, 'EdgeColor','k');hold on;

        % plot(cell_verts(:,1), cell_verts(:,2), 'ro','MarkerSize',5);hold on;
        % Cx = mean(cell_verts(:,1)); Cy = mean(cell_verts(:,2));
        % if (Cx >= 0 && Cx <= bx) && (Cy >= 0 && Cy <= bx)
        %    text(Cx, Cy-1, num2str(cell_no), 'Color','k','FontSize', 6, 'HorizontalAlignment', 'center');
        % end
    end

    % all_cell_verts = unique(all_cell_verts,"rows");
    % plot(all_cell_verts(:,1), all_cell_verts(:,2), 'ro','MarkerSize',5);hold on;

    voronoi_edges = [];

    for k = 1:size(C, 1)
        cell_vertex_indices = C{k};
        cell_verts = V(cell_vertex_indices, :);
        for l = 1:length(cell_vertex_indices)
            if l < length(cell_vertex_indices)
                v1 = cell_verts(l, :);
                v2 = cell_verts(l + 1, :);
            else
                v1 = cell_verts(l, :);
                v2 = cell_verts(1, :);
            end

            [xi, yi] = polyxpoly([v1(1), v2(1)], [v1(2), v2(2)], ...
                constraints(:,1), constraints(:,2));

            if ~isempty(xi)
                if any(v1 < 0) || any(v1 > bx)
                    v1 = [xi(1), yi(1)];
                end
                if any(v2 < 0) || any(v2 > bx)
                    v2 = [xi(1), yi(1)];
                end
                if length(xi) == 2
                    v1 = [xi(1), yi(1)];
                    v2 = [xi(2), yi(2)];
                end
            end

            edge = [v1, v2];
            voronoi_edges = [voronoi_edges; edge];
        end
    end

    edges_info = zeros(size(voronoi_edges, 1), 2);
    labeled_edges = zeros(0, 4);

    for i = 1:size(voronoi_edges, 1)
        x1 = voronoi_edges(i, 1);
        y1 = voronoi_edges(i, 2);
        x2 = voronoi_edges(i, 3);
        y2 = voronoi_edges(i, 4);

        if any(all(labeled_edges == [x1, y1, x2, y2], 2)) || ...
                any(all(labeled_edges == [x2, y2, x1, y1], 2))
            continue
        end

        labeled_edges = [labeled_edges; [x1, y1, x2, y2]];
        length = sqrt((x2 - x1)^2 + (y2 - y1)^2);

        angle = atan2d(y2 - y1, x2 - x1);
        if angle < 0
            angle = angle + 360;
        end

        edges_info(i, :) = [length, angle];
        mid_x = (x1 + x2) / 2;
        mid_y = (y1 + y2) / 2;

        % % Display length and angle
        % text(mid_x, mid_y, sprintf('L=%.2f, A=%.1f°', length, angle), ...
        %     'Color', 'r', 'FontSize', 8, 'FontWeight', 'normal', ...
        %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        % hold on;
    end


    edges_matrix = [voronoi_edges, edges_info];
    edges_matrix = unique(round(edges_matrix, 6), 'rows');
    all_cell_verts = round(all_cell_verts,6);

    D_nodes = unique(all_cell_verts,"rows");

    dist_matrix = pdist2(D_nodes, D_nodes);
    n_points = size(D_nodes, 1);
    conn_matrix = zeros(n_points);

    for i = 1:size(edges_matrix, 1)
        node1 = find(ismember(D_nodes, edges_matrix(i, 1:2), 'rows'));
        node2 = find(ismember(D_nodes, edges_matrix(i, 3:4), 'rows'));
        if ~isempty(node1) && ~isempty(node2)
            conn_matrix(node1, node2) = 1;
            conn_matrix(node2, node1) = 1;
        end
    end

    % gplot(conn_matrix, D_nodes, ':k'); hold on;
    % xlim([0 bx]); ylim([0 bx]);
    % xticks([]); yticks([]);
    % set(gcf, 'Renderer', 'opengl');

    % for i = 1:size(D_nodes, 1)
    %     x = D_nodes(i, 1); y = D_nodes(i, 2);
    %     text(x, y, num2str(i), 'Color', 'blue', 'FontSize', 8, ...
    %          'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    % end

    lowest_angle_matrix = zeros(size(D_nodes, 1), 2);

    for i = 1:size(D_nodes, 1)

        if round(D_nodes(i, 1),4) == bx || round(D_nodes(i, 2),4) == 0 || D_nodes(i, 2) == bx
            continue;
        end

        connected_nodes = find(conn_matrix(i, :) == 1);
        %      plot(D_nodes(i,1),D_nodes(i,2),'rs');hold on;

        if isempty(connected_nodes)
            continue;
        end

        min_angle = Inf;
        min_angle_node = -1;

        valid_angles = [];
        valid_nodes = [];

        for j = 1:size(connected_nodes, 2)
            connected_node = connected_nodes(j);

            if connected_node == i
                continue;
            end

            x1 = D_nodes(i, 1); y1 = D_nodes(i, 2);
            x2 = D_nodes(connected_node, 1); y2 = D_nodes(connected_node, 2);

            if y2 == bx || y2 == 0 || x2 == 0
                continue;
            end

            angle = atan2d((y2 - y1), (x2 - x1));

            if abs(round(angle, 4)) >= 90
                continue;
            end

            valid_angles = [valid_angles; angle];
            valid_nodes = [valid_nodes; connected_node];
        end


        if ~isempty(valid_angles)

            [~, sorted_idx] = sort(abs(valid_angles));
            min_angle = valid_angles(sorted_idx(1));
            min_angle_node = valid_nodes(sorted_idx(1));

            x1 = D_nodes(i, 1); y1 = D_nodes(i, 2);
            x2 = D_nodes(min_angle_node, 1); y2 = D_nodes(min_angle_node, 2);
            % line([x1, x2], [y1, y2], 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '-');
        end

        lowest_angle_matrix(i, :) = [i,abs(round(min_angle,2))];
    end

    lowest_angle_matrix(lowest_angle_matrix(:, 1) == 0, :) = [];
    lowest_angle_matrix(lowest_angle_matrix(:, 2) == Inf, :) = [];

    xlim([0 bx]); ylim([0 bx]);
    xticks([]); yticks([]);
    set(gcf, 'Renderer', 'opengl');

    A = lowest_angle_matrix(:, 2);
    % disp(A);

end
A = lowest_angle_matrix(:, 2);
% figure;
% histogram(A, 50, 'Normalization', 'pdf');hold on;
% xlabel('theta');
% ylabel('PDF');


theta = linspace(min(A), max(A), 500); 
[counts, edges] = histcounts(A, theta, 'Normalization', 'pdf');
centers = edges(1:end-1) + diff(edges)/2; 
focus_region = centers >= 80 & centers <= 90; 
theta_focus = centers(focus_region);
pdf_focus = counts(focus_region);
fit_func = @(params, theta) params(1) * (1 - theta / 90).^params(2);
initial_guess = [0.06, 2];
options = optimoptions('lsqcurvefit', 'Display', 'iter');
lb = [0, 0]; % Lower bounds for C and n
ub = [Inf, Inf]; % Upper bounds for C and n
params = lsqcurvefit(fit_func, initial_guess, theta_focus, pdf_focus, lb, ub, options);

C = params(1);
n = params(2);

figure;
histogram(A, 100, 'Normalization', 'pdf'); hold on;
plot(theta_focus, fit_func(params, theta_focus), 'r-', 'LineWidth', 2);
xlabel('\theta');
ylabel('PDF');
legend('Histogram', sprintf('Fit: C(1-\\theta/90)^n\nC=%.2f, n=%.2f', C, n));
title('PDF Fit Near \theta = 90');
hold off;






