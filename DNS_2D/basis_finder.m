function [K_IC,rotated_planes,planes,sortedIndices,SIF,path,allBasis] = basis_finder(system_rotation_angle)
% if system_rotation_angle > 90
%     system_rotation_angle = 90 - system_rotation_angle;
% end
path = [];

% ________________________________Eight fold symmetry
% hundred_type_planes = [1 0;-1 0;0 1;0 -1];
% one_ten_type_planes = [1 1;1 -1;-1 -1;-1 1;
%              1/2  1;1/2 -1;-1/2 -1;-1/2 1;
%              1/4  1;1/4 -1;-1/4 -1;-1/4 1;
%              1 1/2;1 -1/2;-1 -1/2;-1 1/2;
%              1 1/4;1 -1/4;-1 -1/4;-1 1/4];

% % ________________________________four fold symmetry
hundred_type_planes = [1 0;-1 0;0 1;0 -1];
one_ten_type_planes = [1 1;1 -1;-1 -1;-1 1];

% ________________________________three fold symmetry
% hundred_type_planes = [1 0;-1 0];
% one_ten_type_planes = [sind(60) 1;sind(60) -1;-sind(60) -1;-sind(60) 1];

% ________________________________two fold symmetry
% hundred_type_planes = [1 0;-1 0;0 1;0 -1];
% one_ten_type_planes = [];

% ________________________________one fold symmetry
% hundred_type_planes = [1 0;-1 0];
% one_ten_type_planes = [];


planes = [hundred_type_planes; one_ten_type_planes];
rotation_theta = deg2rad(system_rotation_angle);

R_theta = [cos(rotation_theta), -sin(rotation_theta);
    sin(rotation_theta), cos(rotation_theta)];
rotated_planes = (R_theta * planes')';

theta = zeros(size(planes, 1), 1);
for i = 1:size(rotated_planes, 1)
    vec = rotated_planes(i, :);
    theta(i) = rad2deg(atan2(vec(2), vec(1)));
end

start_index_110 = size(hundred_type_planes, 1) + 1;
end_index_110 = start_index_110 + size(one_ten_type_planes, 1) - 1;

CE_100 = 3.00; CE_110 = 6.48;
CE = zeros(size(planes, 1), 1);
CE(1:start_index_110-1) = CE_100;             % For {100}
CE(start_index_110:end_index_110) = CE_110;   % For {110}

E = 568*10^9; nu = 0.25;

K_IC = zeros(size(planes, 1), 1);
allBasis = cell(size(planes, 1), 1);
color(1:(start_index_110-1))='r'; color(start_index_110:end_index_110)='b';


for i = 1:size(rotated_planes, 1)
    basis1 = rotated_planes(i, :);
    basis1 = basis1/norm(basis1);
    if basis1(1) < 0
        basis1 = -basis1;
    end
    basis0 = planes(i,:);
    basis0 = basis0/norm(basis0);
    allBasis{i} = basis1;
    theta1 = rad2deg(atan2(basis1(2), basis1(1)));
    K_IC(i) = (1) ./ ((cosd(theta1)).^1);
end

[~, sortedIndices] = sort(K_IC, 'ascend');
SIF = K_IC(sortedIndices);
paths = rotated_planes(sortedIndices, :);
% allBasis = allBasis(sortedIndices);

path_basis = paths(1,:)/norm(paths(1,:));
path.KIC = SIF(1);path.current = path_basis;
path.parent = planes(sortedIndices(1),:);
path.angle_x = rad2deg(atan2(path_basis(2), path_basis(1)));

end



