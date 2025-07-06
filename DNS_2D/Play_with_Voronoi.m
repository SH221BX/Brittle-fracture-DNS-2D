clear all;
rng(1,"v4");
bx = 100;
constraints=[0 0; 0 1; 1 1; 1 0; 0 0] * bx;

% num_columns = 7; 
% cell_width = bx/num_columns; margin = cell_width/2;
% centerline = bx / 2;
% x_positions = linspace(margin, bx-margin, num_columns); 
% grid_points = [x_positions', repmat(centerline, num_columns, 1)];


num_rows = 7; 
cell_width = bx/num_rows; margin = cell_width/2;
centerline = bx / 2;
y_positions = linspace(margin, bx-margin, num_rows); 
grid_points = [repmat(centerline, num_rows, 1),y_positions'];


[V, C] = VoronoiBounded(grid_points(:,1), grid_points(:,2), constraints);

figure; hold on; axis equal; box on;
colormap('jet'); 

for cell_no = 1:length(C)
    if all(C{cell_no}~=1) 
        verts = V(C{cell_no},:);
        fill(verts(:,1), verts(:,2), rand(1,3), 'FaceAlpha', 0, 'EdgeColor', 'k'); % Fill cells with random color
    end
end

plot(grid_points(:,1), grid_points(:,2), 'ks', 'MarkerFaceColor', 'w');

xticks([]); yticks([]); 
xlim([0 bx]); ylim([0 bx]); 
set(gcf, 'Renderer', 'opengl'); 


%%
clear all;

bx = 100; 
constraints=[0 0; 0 1; 1 1; 1 0;0 0] * bx;

aspect_ratio = 1; 
rows = 5; 
cols = rows * aspect_ratio; 

cell_width = bx / cols; cell_height = bx / rows;
margin_x = cell_width / 2; margin_y = cell_height / 2;

x = linspace(margin_x, bx - margin_x, cols);
y = linspace(margin_y, bx - margin_y, rows);
[X, Y] = meshgrid(x, y);
grid_points = [X(:), Y(:)];

[V, C] = VoronoiBounded(grid_points(:,1), grid_points(:,2), constraints);

figure;hold on;
axis equal; box on;
colormap('jet'); 

for cell_no = 1:length(C)
    if all(C{cell_no}~=1)
        verts = V(C{cell_no},:);
        fill(verts(:,1), verts(:,2), rand(1,3),'FaceAlpha',0, 'EdgeColor', 'k'); % Fill cells with random color
    end
end

plot(grid_points(:,1), grid_points(:,2), 'ks', 'MarkerFaceColor', 'w');


xticks([]); yticks([]); xlabel(''); ylabel('');
xlim([0 bx]); ylim([0 bx]); set(gcf, 'Renderer', 'opengl');





% clear all;
% 
% rng(15,'v4'); GBC = 1;
%     bx = 100; edges = [0 0 1 1 0;0 1 1 0 0] * bx;
%     num_nodes = 50;
%     num_nodes = num_nodes-4;
%     nodes_pos = rand(num_nodes, 2)*bx;
% 
%     constraints=[0 0; 0 1; 1 1; 1 0;0 0] * bx;
%     nodes_pos = [nodes_pos; constraints]; nodes_pos = unique(nodes_pos, 'rows');
% 
%     x = nodes_pos(:, 1); y = nodes_pos(:, 2);
%     dt = delaunayTriangulation(x,y); [V, C] = voronoiDiagram(dt);
% 
% 
%     cell_orientations = -90 + (180) * rand(1, size(C, 1));
%     all_cell_verts = []; int_points = [];grains_lib = []; colormap('hsv');
% 
%     for cell_no = 1:size(C, 1)
%         [cell_verts] = finding_xo(V, C, cell_no, bx, constraints);
%         scatter(cell_verts(:,1),cell_verts(:,2),30);hold on;
%         centroid = mean(cell_verts);
%         sorting_angles1 = atan2(cell_verts(:,2) - centroid(2), cell_verts(:,1) - centroid(1));
%         [~, order1] = sort(sorting_angles1); 
%         cell_verts = cell_verts(order1, :);
%         all_cell_verts = [all_cell_verts;cell_verts];
%         orientation = cell_orientations(cell_no);
%         color_value = orientation;
%         patch(cell_verts(:,1), cell_verts(:,2),color_value,'FaceAlpha',0, 'EdgeColor','k');hold on;
% 
%         plot(cell_verts(:,1), cell_verts(:,2), '-k');hold on;
%     end
% 
% 
%     xticks([]);yticks([]);xlabel('');ylabel(''); axis equal;
%     grid off;box on; xlim([0 bx]); ylim([0 bx]);
%     %xlabel('X','Interpreter','latex'); ylabel('Y','Interpreter','latex');
%     set(gcf, 'Renderer', 'opengl');


%%

clear all;

bx = 100; 
constraints=[0 0; 0 1; 1 1; 1 0;0 0] * bx;

% Hexagonal Grains________________________________
side_length = 6;
hex_height = 1.5*sqrt(2) * side_length;
hex_width = 5*sqrt(2)* side_length;
num_cols = ceil((130) / (1.5 * side_length)) + 2;
num_rows = ceil((130) / hex_height) + 2;
nodes_pos = [];

for row = 1:num_rows
    for col = 1:num_cols
        x = (col - 1) * 1.5 * side_length - 20; % Start from -20
        if mod(row, 2) == 0
            x = x + 0.75 * side_length;
        end
        y = (row - 1) * hex_height - 20; % Start from -20
        nodes_pos = [nodes_pos; x, y];
    end
end

% Adjusting the boundary [-20, 110]  _ _ _ _ _
margin = 10;
nodes_pos = nodes_pos(nodes_pos(:,1) >= (-20 - margin) & nodes_pos(:,1) <= (110 + margin) & ...
nodes_pos(:,2) >= (-20 - margin) & nodes_pos(:,2) <= (110 + margin), :);

[V, C] = VoronoiBounded(nodes_pos(:,1), nodes_pos(:,2), constraints);

figure;hold on;
axis equal; box on;
colormap('jet'); 

for cell_no = 1:length(C)
    if all(C{cell_no}~=1)
        verts = V(C{cell_no},:);
        fill(verts(:,1), verts(:,2), rand(1,3),'FaceAlpha',0, 'EdgeColor', 'k'); % Fill cells with random color
    end
end

plot(nodes_pos(:,1), nodes_pos(:,2), 'ks', 'MarkerFaceColor', 'w');


xticks([]); yticks([]); xlabel(''); ylabel('');
xlim([0 bx]); ylim([0 bx]); set(gcf, 'Renderer', 'opengl');


%%

clear all;

brick_width = 10;
brick_height = 5;
num_rows = 20; 
num_cols = 10;

% Figure setup
figure;
hold on;

% Loop over rows and columns to draw bricks
for row = 0:num_rows-1
    for col = 0:num_cols-1
        % Calculate the position of the brick
        x_offset = col * brick_width;
        y_offset = row * brick_height;

        % For every second row, offset the bricks by half a brick width
        if mod(row, 2) == 1
            x_offset = x_offset + brick_width / 2;
        end

        % Draw the brick as a rectangle
        rectangle('Position', [x_offset, y_offset, brick_width, brick_height], ...
                  'FaceColor', [1, 0.8, 0.6], 'EdgeColor', [0.8, 0.4, 0.1]);
    end
end

% Axis settings
axis equal; xticks([]); yticks([]); box on;
xlim([0 num_cols * brick_width + brick_width / 2]);
ylim([0 num_rows * brick_height]);
title('Brick and Mortar Pattern');
hold off;
%%
clear all;

x = 10;
brick_width = x;
brick_height = x/2;
bx = 100;  

num_cols = floor(bx / brick_width);
num_rows = floor(bx / brick_height);

figure; hold on; axis equal; box on;


for row = 0:num_rows-1
    y_offset = row * brick_height;
    
    if mod(row, 2) == 0
        x_offset = -brick_width / 2;
        num_bricks_in_row = ceil(bx / brick_width) + 1; 
    else
        x_offset = 0;
        num_bricks_in_row = ceil(bx / brick_width);
    end

    for col = 0:num_bricks_in_row-1
        current_brick_width = brick_width;

        if x_offset < 0  
            current_brick_width = brick_width + x_offset;
            x_offset = 0;  
        end

        if x_offset + current_brick_width > bx  % Adjust for right edge
            current_brick_width = bx - x_offset;  % Adjust to not exceed the boundary
        end

        rectangle('Position', [x_offset, y_offset, current_brick_width, brick_height], ...
                  'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1);

        x_offset = x_offset + current_brick_width;  % Increment x_offset
    end
end

xticks([]); yticks([]);
xlim([0, bx]); 
ylim([0, bx]);  
set(gcf, 'Renderer', 'opengl');
%%


clear all;

bx = 100;
constraints=[0 0; 0 1; 1 1; 1 0; 0 0] * bx;
x = 20;
brick_width = x;
brick_height = x / 2;
bx = 100;

num_cols = floor(bx / brick_width);
num_rows = floor(bx / brick_height);

figure; hold on; axis equal; box on;

V = [];  C = {}; 

vertex_count = 0;

for row = 0:num_rows-1
    y_offset = row * brick_height;
    
    if mod(row, 2) == 0
        x_offset = -brick_width / 2;
        num_bricks_in_row = ceil(bx / brick_width) + 1; 
    else
        x_offset = 0;
        num_bricks_in_row = ceil(bx / brick_width);
    end

    for col = 0:num_bricks_in_row-1
        current_brick_width = brick_width;

        if x_offset < 0  
            current_brick_width = brick_width + x_offset;
            x_offset = 0;  
        end

        if x_offset + current_brick_width > bx  % Adjust for right edge
            current_brick_width = bx - x_offset;  % Adjust to not exceed the boundary
        end

        rectangle('Position', [x_offset, y_offset, current_brick_width, brick_height], ...
                  'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1);

        V = [V; x_offset, y_offset; 
                x_offset + current_brick_width, y_offset; 
                x_offset + current_brick_width, y_offset + brick_height; 
                x_offset, y_offset + brick_height];
        
        % Add the indices plus the first index again at the end to close the loop
        C{end+1} = [vertex_count + 1, vertex_count + 2, vertex_count + 3, vertex_count + 4, vertex_count + 1];
       
        vertex_count = vertex_count + 4;
        x_offset = x_offset + current_brick_width;  
    end
end

xticks([]); yticks([]);
xlim([0, bx]); ylim([0, bx]);  
set(gcf, 'Renderer', 'opengl');

% Display the vertices and closed loops
figure;
for cell_no = 1:length(C)
    cell_verts = V(C{cell_no}(1:end-1), :); % Avoid repeating the last point for scatter
    scatter(cell_verts(:,1), cell_verts(:,2), 50, 'rs'); hold on;
    centroid = mean(cell_verts);
    sorting_angles = atan2(cell_verts(:,2) - centroid(2), cell_verts(:,1) - centroid(1));
    [~, order] = sort(sorting_angles); 
    cell_verts = [cell_verts(order, :); cell_verts(order(1),:)];
    plot(cell_verts(:,1), cell_verts(:,2), '-k'); hold on; % Draw closed loop
end

xticks([]); yticks([]);
xlim([0, bx]); ylim([0, bx]);  
set(gcf, 'Renderer', 'opengl');  



%%

clear all;

bx = 100;
x = 20;
brick_width = x;
brick_height = x / 2;

num_rows = floor(bx / brick_width);  
num_cols = floor(bx / brick_height);  

figure; hold on; axis equal; box on;


V = [];   C = {};  
vertex_count = 0;

for col = 0:num_cols-1
    x_offset = col * brick_height;  % This was originally y_offset

    if mod(col, 2) == 0
        y_offset = -brick_width / 2;  % This was originally x_offset
        num_bricks_in_col = ceil(bx / brick_width) + 1;
    else
        y_offset = 0;
        num_bricks_in_col = ceil(bx / brick_width);
    end

    for row = 0:num_bricks_in_col-1
        current_brick_height = brick_width;  

        if y_offset < 0  
            current_brick_height = brick_width + y_offset;
            y_offset = 0;  
        end

        if y_offset + current_brick_height > bx  
            current_brick_height = bx - y_offset;
        end

        rectangle('Position', [x_offset, y_offset, brick_height, current_brick_height], ...
                  'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1);

        V = [V; x_offset, y_offset; 
                x_offset + brick_height, y_offset; 
                x_offset + brick_height, y_offset + current_brick_height; 
                x_offset, y_offset + current_brick_height];
        
        C{end+1} = vertex_count + (1:4);
        vertex_count = vertex_count + 4;
        
        y_offset = y_offset + current_brick_height;  % Increment y_offset
    end
end

xticks([]); yticks([]);
xlim([0, bx]); 
ylim([0, bx]);
set(gcf, 'Renderer', 'opengl');

%%

clear all;
bx = 100;
theta = pi/8;
x = 10;
brick_width = x;
brick_height = x / 2;

num_cols = floor(bx / brick_width);
num_rows = floor(bx / brick_height);

cx = bx / 2;
cy = bx / 2;

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

figure; hold on; axis equal; box on;

V = [];  C = {}; 
vertex_count = 0;

for row = 0:num_rows-1
    for col = 0:num_cols-1
        y_offset = row * brick_height;
        x_offset = col * brick_width - (mod(row, 2) * brick_width / 2);

        current_brick_width = brick_width;

        if x_offset + brick_width > bx
            current_brick_width = bx - x_offset;
        end

        vertices = [x_offset, y_offset;
                    x_offset + current_brick_width, y_offset;
                    x_offset + current_brick_width, y_offset + brick_height;
                    x_offset, y_offset + brick_height];

        vertices = vertices - [cx, cy]; 
        vertices = (R * vertices')'; 
        vertices = vertices + [cx, cy]; 

        V = [V; vertices];
        C{end + 1} = vertex_count + (1:4);
        vertex_count = vertex_count + 4;

        fill(vertices(:,1), vertices(:,2), 'w', 'EdgeColor', 'k', 'LineWidth', 1);
    end
end

% Display vertices and cells graphically
% figure;
% for cell_no = 1:length(C)
%     cell_verts = V(C{cell_no},:);
%     scatter(cell_verts(:,1), cell_verts(:,2), 50, 'rs'); hold on;
%     centroid = mean(cell_verts);
%     sorting_angles = atan2(cell_verts(:,2) - centroid(2), cell_verts(:,1) - centroid(1));
%     [~, order] = sort(sorting_angles);
%     cell_verts = [cell_verts(order, :); cell_verts(order(1),:)];
% %     plot(cell_verts(:,1), cell_verts(:,2), '-k'); hold on;
% end

xticks([]); yticks([]);
xlim([0, bx]); ylim([0, bx]);
set(gcf, 'Renderer', 'opengl');

%%
clear all;

bx = 100;                % Bounding box size
theta = pi/4;            % Rotation angle (45 degrees)
brick_width = 25;        % Brick width
brick_height = 10;        % Brick height

diagonal = sqrt(2 * (bx ^ 2))+1;
num_cols = ceil(diagonal / brick_width) + 2;  
num_rows = ceil(diagonal / brick_height) + 2; 

cx = bx / 2;
cy = bx / 2;

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
V = [];
C = [];
vertex_count = 0;

figure; hold on; axis equal; box on;


for row = 0:num_rows-1
    for col = 0:num_cols-1
        y_offset = row * brick_height - (num_rows * brick_height / 2 - cy);
        x_offset = col * brick_width - (num_cols * brick_width / 2 - cx) - (mod(row, 2) * brick_width / 2);

        vertices = [x_offset, y_offset;
                    x_offset + brick_width, y_offset;
                    x_offset + brick_width, y_offset + brick_height;
                    x_offset, y_offset + brick_height];

   
        vertices = (R * (vertices - [cx, cy])')' + [cx, cy];

        V = [V; vertices];
        C = [C; vertex_count + (1:4)];
        vertex_count = vertex_count + 4;

        if any(inpolygon(vertices(:,1), vertices(:,2), [0 bx bx 0], [0 0 bx bx]))
            patch('Vertices', vertices, 'Faces', [1 2 3 4], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
end

plot([0 bx bx 0 0], [0 0 bx bx 0], 'r-', 'LineWidth', 2);

xticks([]); yticks([]);
xlim([0, bx]);
ylim([0, bx]);
set(gcf, 'Renderer', 'opengl');



%%


clear all;

bx = 100;                % Bounding box size
theta = pi/4;            % Rotation angle (45 degrees)
brick_width = 10;        % Brick width
brick_height = 5;        % Brick height

% Compute grid size for full coverage
diagonal = sqrt(2 * (bx ^ 2));
num_cols = ceil(diagonal / brick_width) + 2;  % Add two extra columns
num_rows = ceil(diagonal / brick_height) + 2;  % Add two extra rows

cx = bx / 2;  % Center of rotation
cy = bx / 2;

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];  % Rotation matrix

figure; hold on; axis equal; box on;
title('Rotated Brick Pattern at 45 Degrees');

for row = 0:num_rows-1
    for col = 0:num_cols-1
        x_offset = col * brick_width - (num_cols * brick_width / 2) + cx;
        y_offset = row * brick_height - (num_rows * brick_height / 2) + cy;

        % Define vertices
        vertices = [x_offset, y_offset;
                    x_offset + brick_width, y_offset;
                    x_offset + brick_width, y_offset + brick_height;
                    x_offset, y_offset + brick_height];

        % Rotate vertices
        vertices = (R * (vertices - [cx, cy])')' + [cx, cy];

        % Clipping each vertex
        for i = 1:size(vertices, 1)
            if vertices(i, 1) < 0
                vertices(i, 1) = 0;
            elseif vertices(i, 1) > bx
                vertices(i, 1) = bx;
            end
            if vertices(i, 2) < 0
                vertices(i, 2) = 0;
            elseif vertices(i, 2) > bx
                vertices(i, 2) = bx;
            end
        end

        % Draw the brick
        fill(vertices(:,1), vertices(:,2), 'w', 'EdgeColor', 'k', 'LineWidth', 1);
    end
end

% Draw the bounding box
plot([0 bx bx 0 0], [0 0 bx bx 0], 'r-', 'LineWidth', 2);

xticks([]); yticks([]);
xlim([0, bx]);
ylim([0, bx]);
set(gcf, 'Renderer', 'opengl');

