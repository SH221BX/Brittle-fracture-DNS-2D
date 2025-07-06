
clear all; clc; format long;

R0 = 1; R = 1;
n_segments = 1000;

A = 1.09530; M = -12.2242; phi = 0.8483;   
R_max = R/cos(pi/8);
z_min = 1.00001; z_max = R_max;
X = 1000;
z = linspace(z_min, z_max, X);

f_KV = @(x) (A ./ (1 + exp(M * ((pi/2 - acos(R0 ./x)).^(1/3) - phi))) ) .* (R0 ./ x.^2) ./ sqrt(1 - (R0 ./x).^2) .* (x >= R & x <= R / cos(pi/8));
f_KT = @(x) (8 * R) ./ (pi * x.^2 .* sqrt(1 - (R ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));
f_KG = @(x) (2 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));

F_KV = arrayfun(@(xx) integral(f_KV, R0, xx), z);
F_KT = arrayfun(@(xx) integral(f_KT, R, min(xx, R ./ cos(pi/8))), z);
F_KG = arrayfun(@(xx) integral(f_KG, R0, xx), z);

fZ1 = f_KV(z) .* (1 - F_KT) + f_KT(z) .* (1 - F_KV);
fZ2 = f_KG(z) .* (1 - F_KT) + f_KT(z) .* (1 - F_KG);

F_Z1 = cumtrapz(z, fZ1); F_Z1 = F_Z1 / max(F_Z1);
F_Z2 = cumtrapz(z, fZ2); F_Z2 = F_Z2 / max(F_Z2);

f_S = cell(1, n_segments); F_S = cell(1, n_segments);
F_Max = cell(1, n_segments); f_Max = cell(1, n_segments);

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
    fK = f_Max{seg};
end

% filename = 'R1_two_draw_vor.gif';
% fig = figure('Color', 'w');
% % set(fig, 'Position', [200, 300, 400, 300]);
% 
% for seg = 1:n_segments    
%     subplot(2,ceil(n_segments/2),seg);
%     plot(z,f_Max{seg},'r-','LineWidth',1.5);
%     title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
%     set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
% 
%     grid off; hold off;
%     drawnow; pause(0);
% 
%     % frame = getframe(gcf);
%     % img = frame2im(frame);
%     % [imind, cm] = rgb2ind(img, 256);
%     % if seg == 1 
%     %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
%     % else 
%     %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
%     % end
% 
% end


F_K = cumtrapz(z, fK);
F_K = F_K / max(F_K);
Qvals = 1:1:100;

% figure;
colors = lines(length(Qvals));

for i = 1:1:100
    Q = Qvals(i);
    % subplot(2, ceil(length(Qvals)/2), i);
    f_M = Q .* (F_K).^(Q-1) .* fK;
    % plot(z, f_M, 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('Q = %d', Q));
    % title(['$\rm Q = \, ' num2str(Q) '$'],'Interpreter','latex', 'FontSize',8);
    % set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % grid off;

    norm_factor = trapz(z,f_M);
    mean_pdf(i) = trapz(z,z.*f_M)/norm_factor;
end

% figure;
plot(1:100,mean_pdf,'-','LineWidth',2); hold on;
xlabel('N','Interpreter','latex','FontSize',14);
ylabel('$\mathrm{E[K_{M}]}$','Interpreter','latex','FontSize',14);
set(gca, 'FontSize', 14, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
ylim([1.063 1.0825]);



%%


clear; clc; close all;

N = 1000000; 
% theta = acos(1 - 2*rand(N,1)); 
theta   = 2*pi*rand(N,1);

phi   = 2*pi*rand(N,1);

% 1) Distributions of theta and phi
% Note: For uniform normals on the sphere:
%   theta in [0, pi],  pdf(theta) = (1/2)*sin(theta)
%   phi   in [0,2*pi], pdf(phi)   = 1/(2*pi)

figure('Color','w');
subplot(1,2,1);
histogram(theta, 'Normalization','pdf','BinEdges', linspace(0, pi, 50),...
    'FaceColor',[0.2 0.6 0.8],'EdgeColor','w');
hold on;
tt = linspace(0,pi,200);
% plot(tt, 0.5*sin(tt),'r-','LineWidth',2,'DisplayName','$\frac12 \sin(\theta)$');
xlabel('$\theta$','Interpreter','latex','FontSize',14);
ylabel('PDF','Interpreter','latex','FontSize',14);
title('Distribution of $\theta$','Interpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','north');
grid on;

subplot(1,2,2);
histogram(phi, 'Normalization','pdf','BinEdges', linspace(0,2*pi,50),...
    'FaceColor',[0.9 0.4 0.5],'EdgeColor','w');
hold on;
pp = linspace(0,2*pi,200);
% plot(pp, 1/(2*pi)*ones(size(pp)),'r-','LineWidth',2,'DisplayName','$\frac{1}{2\pi}$');
xlabel('$\phi$','Interpreter','latex','FontSize',14);
ylabel('PDF','Interpreter','latex','FontSize',14);
title('Distribution of $\phi$','Interpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','north');
grid on;

N_subset = 20;
indices = randperm(N, N_subset);

figure('Color','w');
[xs,ys,zs] = sphere(40);
surf(xs,ys,zs,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',[0.2 0.8 0.8]);
hold on;
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');

for i = 1:N_subset
    th = theta(indices(i));
    ph = phi(indices(i));
    n = [sin(th)*cos(ph); sin(th)*sin(ph); cos(th)];

    if abs(n(1))<0.9
        temp = [1;0;0];
    else
        temp = [0;1;0];
    end
    e2 = cross(n,temp); e2 = e2/norm(e2);
    e3 = cross(n,e2);   e3 = e3/norm(e3);
    
    svec = linspace(0,2*pi,100);
    circle = zeros(3,numel(svec));
    for k = 1:numel(svec)
        circle(:,k) = e2*cos(svec(k)) + e3*sin(svec(k));
    end
    plot3(circle(1,:), circle(2,:), circle(3,:), 'LineWidth',1);
    drawnow;
end
hold off;


%%



clear; clc; close all;

N = 1e7; 
xyz = randn(N,3);
r = sqrt(sum(xyz.^2,2));
xyz_unit = xyz ./ r;

x = xyz_unit(:,1);
y = xyz_unit(:,2);
z = xyz_unit(:,3);

theta = acos(z);
phi   = atan2(y,x);
phi(phi<0) = phi(phi<0) + 2*pi; 

figure('Color','w');
subplot(1,2,1);
histogram(theta,'Normalization','pdf','BinEdges',linspace(0,pi,60),...
    'FaceColor',[0.2,0.6,0.8],'EdgeColor','w');
hold on;
th_ = linspace(0,pi,200);
pdf_th = 0.5*sin(th_);
plot(th_, pdf_th, 'r--','LineWidth',2,'DisplayName','$\frac12 \sin(\theta)$');
xlabel('$\theta$','Interpreter','latex','FontSize',14);
ylabel('PDF','Interpreter','latex','FontSize',14);
title('Distribution of $\theta$','Interpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','north');
grid on;

subplot(1,2,2);
histogram(phi,'Normalization','pdf','BinEdges',linspace(0,2*pi,60),...
    'FaceColor',[0.9,0.4,0.5],'EdgeColor','w');
hold on;
ph_ = linspace(0,2*pi,200);
pdf_ph = 1/(2*pi)*ones(size(ph_));
plot(ph_, pdf_ph, 'r--','LineWidth',2,'DisplayName','$\frac{1}{2\pi}$');
xlabel('$\phi$','Interpreter','latex','FontSize',14);
ylabel('PDF','Interpreter','latex','FontSize',14);
title('Distribution of $\phi$','Interpreter','latex','FontSize',16);
legend('Interpreter','latex','Location','north');
grid on;






%%


clear all; clc; close all;
N = 1e5;
K0 = 1;

n_xy0 = [0; 0; 1];
n_yz0 = [1; 0; 0];
n_zx0 = [0; 1; 0];

minKeq_vals = zeros(N,1);
z_hat = [0; 0; 1];

for i = 1:N
alpha = 2 * pi * rand;  % Rotation about z, uniform in [0, 2pi]

beta_min = 0;
beta_max = pi/2 ;

cos_beta_min = cos(beta_max);  
cos_beta_max = cos(beta_min);  

cos_beta = cos_beta_min + (cos_beta_max - cos_beta_min) * rand; 
beta = acos(cos_beta); 
gamma = 2 * pi * rand;  
  
    Rz1 = [cos(alpha), -sin(alpha), 0;
           sin(alpha),  cos(alpha), 0;
           0,           0,          1];
       
    Ry = [cos(beta), 0, sin(beta);
          0,         1, 0;
         -sin(beta), 0, cos(beta)];
     
    Rz2 = [cos(gamma), -sin(gamma), 0;
           sin(gamma),  cos(gamma), 0;
           0,           0,          1];
       
    
    R = Rz1 * Ry * Rz2;
   
    n_xy = R * n_xy0;
    n_yz = R * n_yz0;
    n_zx = R * n_zx0;
   
    c_xy = dot(n_xy, z_hat);
    c_yz = dot(n_yz, z_hat);
    c_zx = dot(n_zx, z_hat);
    
    K_eq_xy = K0 / max(eps, abs(c_xy));
    K_eq_yz = K0 / max(eps, abs(c_yz));
    K_eq_zx = K0 / max(eps, abs(c_zx));
    
    K_xy(i) = K_eq_xy;  K_yz(i) = K_eq_yz; K_zx(i) = K_eq_zx;

    minKeq_vals(i) = min([K_eq_xy, K_eq_yz, K_eq_zx]);
end

 % 
 % figure;
 %  histogram(minKeq_vals, 'Normalization','pdf', 'BinEdges', linspace(1,max(minKeq_vals),1000));
 % 

  % Plot the histogram
% figure('Color','w');
% histogram(minKeq_vals, 'Normalization','pdf', 'BinEdges', linspace(1,2,100));
% histogram(K_xy, 'Normalization','pdf', 'BinEdges', linspace(1,2,100));
% figure;
% histogram(K_yz, 'Normalization','pdf', 'BinEdges', linspace(1,2,100));
% figure;
% histogram(K_zx, 'Normalization','pdf', 'BinEdges', linspace(1,2,100));

 
% xlabel('K','FontSize',12);
% ylabel('PDF','FontSize',12);
% % title('Histogram of min(K_{eq}) across xy, yz, zx planes (Euler Angles)','FontSize',14);
% grid on;


figure;

X = acos(1./minKeq_vals);
x = linspace(0, pi/4, 1000);
y = 3*sin(x) ;
plot(x,y,'-k','LineWidth',3);hold on;
x2 = linspace(pi/4, acos(sqrt(1/3)), 1000);
y_start = 3*sin(pi/4);    y_end = 0;               
y2 = y_start + (y_end - y_start) * ((x2 - min(x2)) / (max(x2) - min(x2))).^0.5;
 % y2 = sqrt(2).*(atan(1./sqrt(y2) - pi/3));

histogram(X,'Normalization','pdf','BinEdges',linspace(0,pi/3.2882,200),'FaceColor','c','FaceAlpha',1,'EdgeColor','c');hold on;
plot(x,y,'-k','LineWidth',3);hold on;
plot(x2,y2,'-k','LineWidth',3);hold on;
set(gca, 'FontSize', 20, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

xlim([0 pi/3]);
x_ticks = linspace(0,pi/3, 5);
xticks(x_ticks);
x_tick_labels = {'$0$', '$\frac{\pi}{12}$', '$\frac{\pi}{6}$', '$\frac{\pi}{4}$', '$\frac{\pi}{3}$'};
xticklabels(x_tick_labels);

xlabel('$\rm {\theta_{min}}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('${f_\theta}$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

grid off;
box on;



%%

clear; clc; close all;
N = 1e6;
K0 = 1;

n_xy0 = [0; 0; 1];
n_yz0 = [1; 0; 0];
n_zx0 = [0; 1; 0];

minKeq_vals = zeros(N,1);
K_xy = zeros(N,1); K_yz = zeros(N,1); K_zx = zeros(N,1);
z_hat = [0; 0; 1];

for i = 1:N
    alpha = 2*pi*rand;         % rotation about z, uniform in [0, 2pi]
    beta  = acos(2*rand - 1);    % rotation about y; ensures cos(beta) is uniform in [-1,1]
    gamma = 2*pi*rand;         % second rotation about z, uniform in [0, 2pi]
    
    Rz1 = [cos(alpha), -sin(alpha), 0;
           sin(alpha),  cos(alpha), 0;
           0,           0,          1];
       
    Ry = [cos(beta), 0, sin(beta);
          0,         1, 0;
         -sin(beta), 0, cos(beta)];
     
    Rz2 = [cos(gamma), -sin(gamma), 0;
           sin(gamma),  cos(gamma), 0;
           0,           0,          1];
    
    R = Rz1 * Ry * Rz2;
    
    n_xy = R * n_xy0;
    n_yz = R * n_yz0;
    n_zx = R * n_zx0;
    
    c_xy = dot(n_xy, z_hat);
    c_yz = dot(n_yz, z_hat);
    c_zx = dot(n_zx, z_hat);
    
    K_eq_xy = K0 / max(eps, abs(c_xy));
    K_eq_yz = K0 / max(eps, abs(c_yz));
    K_eq_zx = K0 / max(eps, abs(c_zx));
    
    K_xy(i) = K_eq_xy;
    K_yz(i) = K_eq_yz;
    K_zx(i) = K_eq_zx;
    
    minKeq_vals(i) = min([K_eq_xy, K_eq_yz, K_eq_zx]);
    
end

figure('Color','w'); hold on;
edges = linspace(1,2,100);
histogram(K_xy, 'Normalization','pdf', 'BinEdges', edges, 'FaceColor', 'r', 'FaceAlpha', 0.5);
histogram(K_yz, 'Normalization','pdf', 'BinEdges', edges, 'FaceColor', 'g', 'FaceAlpha', 0.5);
histogram(K_zx, 'Normalization','pdf', 'BinEdges', edges, 'FaceColor', 'b', 'FaceAlpha', 0.5);
histogram(minKeq_vals, 'Normalization','pdf', 'BinEdges', edges, 'FaceColor', 'k', 'FaceAlpha', 0.5);
legend('K_{xy}','K_{yz}','K_{zx}','min(K_{eq})','Location','best');
xlabel('K_{eq}');
ylabel('Probability Density');
title('Histogram of K values for each plane and their minimum');
grid on;





%%

clear all;
phi_min = 0;
phi_max = pi/3.288; 

num_points = 200;
phi_values = linspace(phi_min, phi_max, num_points); 

theta_max = @(phi) acos(1 ./ sqrt((tan(phi).^2 / 2) + 2));
theta_grid = []; phi_grid = [];

for p = phi_values
    max_theta = theta_max(mod(p, pi/2)); 
    theta_values = linspace(0, max_theta, num_points);
    theta_grid = [theta_grid, theta_values];
    phi_grid = [phi_grid, repmat(p, 1, num_points)];
end

% Stereographic Projection Mapping
x = sin(theta_grid) .* cos(phi_grid) ./ (1 + cos(theta_grid));
y = sin(theta_grid) .* sin(phi_grid) ./ (1 + cos(theta_grid));

K = 1 ./ cos(theta_grid);
tri = delaunay(x(:), y(:));
figure;

colormap jet;
trisurf(tri, x(:), y(:), K(:), 'EdgeColor', 'none');
view(2);
axis equal;
colorbar;
set(gca, 'XTick', [], 'YTick', []);
axis off;


plot(x2,y2,'-k','LineWidth',3);hold on;