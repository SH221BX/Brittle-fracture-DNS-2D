clear; clc;

cd('D:\Data_Storage\Toughening_statistics_2D\Com_N1000_1V4\');
Theta_mat1 = load('Theta_mat1.mat', 'Theta_mat').Theta_mat;
Theta_mat2 = load('Theta_mat2.mat', 'Theta_mat').Theta_mat;
KR_mat1 = load('KR_mat1.mat', 'KR_mat').KR_mat;
KR_mat2 = load('KR_mat2.mat', 'KR_mat').KR_mat;

cd('D:\Data_Storage\Toughening_statistics_2D\Com_R12_N1000_1V4\');
KR_mat3 = load('DKR_mat2.mat', 'KR_mat').KR_mat;
cd('C:\Users\sajja\Documents\MATLAB\ASUS\Toughening_statistics_2D');

max_rows = max([size(Theta_mat1, 1), size(Theta_mat2, 1)]);
Theta_mat1 = [Theta_mat1; zeros(max_rows - size(Theta_mat1, 1), size(Theta_mat1, 2))];
Theta_mat2 = [Theta_mat2; zeros(max_rows - size(Theta_mat2, 1), size(Theta_mat2, 2))];
KR_mat1 = [KR_mat1; zeros(max_rows - size(KR_mat1, 1), size(KR_mat1, 2))];
KR_mat2 = [KR_mat2; zeros(max_rows - size(KR_mat2, 1), size(KR_mat2, 2))];
KR_mat3 = [KR_mat3; zeros(max_rows - size(KR_mat3, 1), size(KR_mat3, 2))];


Theta_mat_combined = [Theta_mat1, Theta_mat2];
KR_mat_combined = [KR_mat1, KR_mat2, KR_mat3];
Theta_mat_subset = Theta_mat_combined(1:30, :);
Theta_mat_subset = Theta_mat_subset(~isnan(Theta_mat_subset));
Theta_values = deg2rad(Theta_mat_subset(:));


format long;
N = 1e6;
R0 = 1;
R = 2;
n_segments = 30;

pV = zeros(1, n_segments);
P_T = zeros(1, n_segments);
Segments = cell(1, n_segments);
MaxSegments = cell(1, n_segments);
thetaV = linspace(0, pi/2, N);


 A = 1.052;  M = -12.09;  phi = 0.825;

% fV_unnormalized = 1 ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
% C = trapz(thetaV, fV_unnormalized);
% A = 1 / C;


X = 10000;
z_min = 1.01;
z_max = R/cos(pi/8);
z = linspace(z_min,z_max,X);

f_KV = @(x) (A./(1+exp(M*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).*(x>=R0 & x<R);
f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).*(x>R & x<=R/cos(pi/8));
f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<R);

F_KV = arrayfun(@(xx) integral(f_KV,R0,R),z);
F_KT = arrayfun(@(xx) integral(f_KT,R,min(xx,R/cos(pi/8))),z);
F_KG = arrayfun(@(xx) integral(f_KG,R0,R),z);

F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(z) f_KT(z), R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, R/cos(pi/8));
lower1 = R0; upper1 = R; lower2 = R; upper2 = R/cos(pi/8);
pG_part1 = integral(f_KG, lower1, upper1);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

M = [pV_numerical   pG_numerical
    1-pV_numerical  1-pG_numerical];

P0 = [1; 0];        
           
fZ1 = f_KV(z).*(1-F_KT) + f_KT(z).*(1-F_KV);
fZ2 = f_KG(z).*(1-F_KT) + f_KT(z).*(1-F_KG);
F_Z1 = cumtrapz(z,fZ1); 
F_Z1 = F_Z1/max(F_Z1);
F_Z2 = cumtrapz(z,fZ2); 
F_Z2 = F_Z2/max(F_Z2);

f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
F_Max = cell(1,n_segments);
f_Max = cell(1,n_segments);


for seg = 1:n_segments
    if seg == 1
        P_N = P0;
    else
        [V, D] = eig(M);
        D(1,1) = D(1,1)^(seg-1);
        D(2,2) = D(2,2)^(seg-1);
        M_pow = V * D / V;
        P_N = M_pow * P0;
        disp(P_N);
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

    if seg==1
        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end
end

filename = 'R12_vor_DNS.gif';
fig = figure('Color', 'w','MenuBar','none');
set(fig, 'Position', [200, 300, 400, 300]);

bin_edges = linspace(1, R  / cos(pi/8), 1000);

for seg = 1:n_segments
    % subplot(6, ceil(n_segments/6), seg);
    N = seg;
    KR_sub = KR_mat_combined(1:N, :);
    KR_max = max(KR_sub, [], 1);

    histogram(KR_max, bin_edges, 'Normalization', 'pdf', ...
        'FaceColor', [0, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.8); hold on;

    plot(z, f_Max{seg}, 'r-', 'LineWidth', 2);
    title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    set(gca, 'FontSize', 10, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    ylim([0 8]);
    grid off; hold off;
    drawnow; pause(0.5);

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if seg == 1 
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else 
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end

end



%%


clear; clc;

cd('D:\Data_Storage\Toughening_statistics_2D\Com_R1_N1000_1V4\');
Theta_mat1 = load('Theta_mat1.mat', 'Theta_mat').Theta_mat;
Theta_mat2 = load('Theta_mat2.mat', 'Theta_mat').Theta_mat;
KR_mat1 = load('KR_mat1.mat', 'KR_mat').KR_mat;
KR_mat2 = load('KR_mat2.mat', 'KR_mat').KR_mat;

max_rows = max([size(Theta_mat1, 1), size(Theta_mat2, 1)]);
Theta_mat1 = [Theta_mat1; zeros(max_rows - size(Theta_mat1, 1), size(Theta_mat1, 2))];
Theta_mat2 = [Theta_mat2; zeros(max_rows - size(Theta_mat2, 1), size(Theta_mat2, 2))];
KR_mat1 = [KR_mat1; zeros(max_rows - size(KR_mat1, 1), size(KR_mat1, 2))];
KR_mat2 = [KR_mat2; zeros(max_rows - size(KR_mat2, 1), size(KR_mat2, 2))];

Theta_mat_combined = [Theta_mat1, Theta_mat2];
KR_mat_combined = [KR_mat1, KR_mat2];
Theta_mat_subset = Theta_mat_combined(1:30, :);
Theta_mat_subset = Theta_mat_subset(~isnan(Theta_mat_subset));
Theta_values = deg2rad(Theta_mat_subset(:));

cd('C:\Users\sajja\Documents\MATLAB\ASUS\Toughening_statistics_2D');

R0 = 1; R = 1;
n_segments = 30;

pV = zeros(1, n_segments);
P_T = zeros(1, n_segments);
Segments = cell(1, n_segments);
MaxSegments = cell(1, n_segments);

A = 1.0530; M = -12.2242; phi = 0.8483;

z_min = 1.00001;
z_max = R/cos(pi/8);
R_max = z_max;
X = 1000;
z = linspace(z_min, z_max, X);

f_KV = @(x) (A ./ (1 + exp(M * ((pi/2 - acos(R0 ./x)).^(1/3) - phi))) ) .* (R0 ./ x.^2) ./ sqrt(1 - (R0 ./x).^2) .* (x >= R & x < R / cos(pi/8));
f_KT = @(x) (8 * R) ./ (pi * x.^2 .* sqrt(1 - (R ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));
f_KG = @(x) (2 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (x >= R & x < R / cos(pi/8));

F_KV = arrayfun(@(xx) integral(f_KV, R0, xx), z);
F_KT = arrayfun(@(xx) integral(f_KT, R, min(xx, R ./ cos(pi/8))), z);
F_KG = arrayfun(@(xx) integral(f_KG, R0, xx), z);

fZ1 = f_KV(z) .* (1 - F_KT) + f_KT(z) .* (1 - F_KV);
fZ2 = f_KG(z) .* (1 - F_KT) + f_KT(z) .* (1 - F_KG);
F_Z1 = cumtrapz(z, fZ1);
F_Z1 = F_Z1 / max(F_Z1);
F_Z2 = cumtrapz(z, fZ2);
F_Z2 = F_Z2 / max(F_Z2);
f_S = cell(1, n_segments);
F_S = cell(1, n_segments);
F_Max = cell(1, n_segments);
f_Max = cell(1, n_segments);

F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(z) f_KT(z), R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, R/cos(pi/8));

lower1 = R0; upper1 = R_max; lower2 = R; upper2 = R/cos(pi/8);
pG_part1 = integral(f_KG, lower1, upper1);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

M = [pV_numerical   pG_numerical
    1-pV_numerical  1-pG_numerical];
P0 = [1; 0];  

for seg = 1:n_segments
    if seg == 1
        P_N = P0;
    else
        [V, D] = eig(M);
        D(1,1) = D(1,1)^(seg-1);
        D(2,2) = D(2,2)^(seg-1);
        M_pow = V * D / V;
        P_N = M_pow * P0;
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

    if seg==1
        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end
end

filename = 'R1_vor_DNS.gif';
fig = figure('Color', 'w','MenuBar','none');
set(fig, 'Position', [200, 300, 400, 300]);

bin_edges = linspace(1, R  / cos(pi/8), 100);

for seg = 1:n_segments
    % subplot(6, ceil(n_segments/6), seg);
    N = seg;
    KR_sub = KR_mat_combined(1:N, :);
    KR_max = max(KR_sub, [], 1);

    histogram(KR_max, bin_edges, 'Normalization', 'pdf', ...
        'FaceColor', [0, 0.5, 0.5], 'EdgeColor', 'w', 'FaceAlpha', 0.8); hold on;

    plot(z, f_Max{seg}, 'r-', 'LineWidth', 2);
    title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    set(gca, 'FontSize', 10, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

    grid off; hold off;
    drawnow; pause(0.5);

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if seg == 1 
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else 
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end

end


