%% BASIC of Random Picking -- Empirical K_V and K_G from AK and inclination angles

clear; clc; close all;
rng(1);
N = 1e6;

KIC_GB_raw = readmatrix('AK_data.csv');
theta_IV_raw = readmatrix('angle_selected_rad.csv');
theta_G_raw = readmatrix('angle_all_rad.csv');

KIC_GB_raw = KIC_GB_raw(:);
theta_IV_raw = theta_IV_raw(:);
theta_G_raw = theta_G_raw(:);

KIC_GB_raw = KIC_GB_raw(~isnan(KIC_GB_raw));
theta_IV_raw = theta_IV_raw(~isnan(theta_IV_raw));
theta_G_raw = theta_G_raw(~isnan(theta_G_raw));

KIC_GB_raw = KIC_GB_raw(KIC_GB_raw > 0);
theta_IV_raw = theta_IV_raw(theta_IV_raw >= 0 & theta_IV_raw < pi/2);
theta_G_raw = theta_G_raw(theta_G_raw >= 0 & theta_G_raw < pi/2);

KIC_GB_samples_V = KIC_GB_raw(randi(numel(KIC_GB_raw), N, 1));
theta_IV_samples = theta_IV_raw(randi(numel(theta_IV_raw), N, 1));

K_V_samples = KIC_GB_samples_V ./ cos(theta_IV_samples);

K_V_samples = K_V_samples(isfinite(K_V_samples));
K_V_samples = K_V_samples(K_V_samples > 0);

KIC_GB_samples_G = KIC_GB_raw(randi(numel(KIC_GB_raw), N, 1));
theta_G_samples = theta_G_raw(randi(numel(theta_G_raw), N, 1));

K_G_samples = KIC_GB_samples_G ./ cos(theta_G_samples);

K_G_samples = K_G_samples(isfinite(K_G_samples));
K_G_samples = K_G_samples(K_G_samples > 0);

binEdges_K = 0:1:20;
binEdges_theta = linspace(0, pi/2, 30);

K_V_plot = K_V_samples(K_V_samples >= binEdges_K(1) & K_V_samples <= binEdges_K(end));
K_G_plot = K_G_samples(K_G_samples >= binEdges_K(1) & K_G_samples <= binEdges_K(end));

figure('Color','w','Units','inches','Position',[0.5 0.5 12 6]);

subplot(2,3,1);
hold on; box on; grid off;
histogram(KIC_GB_raw, 40, 'Normalization','pdf', ...
    'FaceColor',[0.2 0.6 1.0], ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$K_{\rm IC}^{\rm GB}$','Interpreter','latex','FontSize',16);
ylabel('$f_{K_{\rm IC}^{\rm GB}}$','Interpreter','latex','FontSize',16);
%title('$K_{\rm IC}^{\rm GB}$','Interpreter','latex','FontSize',15,'FontWeight','normal');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

subplot(2,3,2);
hold on; box on; grid off;
histogram(theta_IV_raw, binEdges_theta, 'Normalization','pdf', ...
    'FaceColor',[0.1 0.35 0.85], ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$\Theta_{\rm I_V}$','Interpreter','latex','FontSize',16);
ylabel('$f_{\Theta}^{\rm I_V}$','Interpreter','latex','FontSize',16);
%title('$\Theta_{\rm I_V}$','Interpreter','latex','FontSize',15,'FontWeight','normal');
xlim([0 pi/2]);
xticks([0 pi/8 pi/4 3*pi/8 pi/2]);
xticklabels({'$0$','$\pi/8$','$\pi/4$','$3\pi/8$','$\pi/2$'});
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

subplot(2,3,3);
hold on; box on; grid off;
histogram(K_V_plot, binEdges_K, 'Normalization','pdf', ...
    'FaceColor',[0.1 0.35 0.85], ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$\rm{K}$','Interpreter','latex','FontSize',16);
ylabel('$f^{\rm{I_V}}_{\rm{K}}$','Interpreter','latex','FontSize',16);
%title('$K_{\rm V}$','Interpreter','latex','FontSize',15,'FontWeight','normal');
xlim([0 20]);
set(gca,'XTick',0:2:20);
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

subplot(2,3,4);
hold on; box on; grid off;
histogram(KIC_GB_raw, 40, 'Normalization','pdf', ...
    'FaceColor',[0.2 0.6 1.0], ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$K_{\rm IC}^{\rm GB}$','Interpreter','latex','FontSize',16);
ylabel('$f_{K_{\rm IC}^{\rm GB}}$','Interpreter','latex','FontSize',16);
%title('$K_{\rm IC}^{\rm GB}$','Interpreter','latex','FontSize',15,'FontWeight','normal');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

subplot(2,3,5);
hold on; box on; grid off;
histogram(theta_G_raw, binEdges_theta, 'Normalization','pdf', ...
    'FaceColor','r', ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$\Theta_{\rm I}$','Interpreter','latex','FontSize',16);
ylabel('$f_{\Theta}^{\rm I}$','Interpreter','latex','FontSize',16);
%title('$\Theta_{\rm I}$','Interpreter','latex','FontSize',15,'FontWeight','normal');
xlim([0 pi/2]);
xticks([0 pi/8 pi/4 3*pi/8 pi/2]);
xticklabels({'$0$','$\pi/8$','$\pi/4$','$3\pi/8$','$\pi/2$'});
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

subplot(2,3,6);
hold on; box on; grid off;
histogram(K_G_plot, binEdges_K, 'Normalization','pdf', ...
    'FaceColor','r', ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$\rm{K}$','Interpreter','latex','FontSize',16);
ylabel('$f^{\rm{I}}_{\rm{K}}$','Interpreter','latex','FontSize',16);
%title('$K_{\rm G}$','Interpreter','latex','FontSize',15,'FontWeight','normal');
xlim([0 20]);
set(gca,'XTick',0:2:20);
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');


%% Finding K_T : EBSD Texture

clear; clc;

K100 = 1.78e6;
K110 = 1.78e6;

planes1 = [1 0 0; 0 1 0; 0 0 1]';
planes2 = [0 1 -1; 0 1 1; 1 -1 0; 1 0 -1; 1 0 1; 1 1 0]';
planes2 = planes2 ./ vecnorm(planes2);

X = load('EBSD_grain_orientation_X3d.mat');
X_3d = X.X_3d;

N = size(X_3d,1);

C1_all = zeros(N,size(planes1,2));
C2_all = zeros(N,size(planes2,2));

K1_vals = zeros(N,1);
K2_vals = zeros(N,1);
minKeq_vals = zeros(N,1);

best_family = strings(N,1);
best_plane_index = zeros(N,1);

theta100_deg = zeros(N,size(planes1,2));
theta110_deg = zeros(N,size(planes2,2));

load_axis = [0 0 1];

for i = 1:N

    phi1 = X_3d(i,1);
    Phi  = X_3d(i,2);
    phi2 = X_3d(i,3);

    Rz1 = [cos(phi1) -sin(phi1) 0; sin(phi1) cos(phi1) 0; 0 0 1];
    Rx  = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
    Rz2 = [cos(phi2) -sin(phi2) 0; sin(phi2) cos(phi2) 0; 0 0 1];

    R = Rz1 * Rx * Rz2;

    planes1_s = R * planes1;
    planes2_s = R * planes2;

    C1 = load_axis * planes1_s;
    C2 = load_axis * planes2_s;

    C1 = abs(C1);
    C2 = abs(C2);

    C1_all(i,:) = C1;
    C2_all(i,:) = C2;

    theta100_deg(i,:) = acosd(C1);
    theta110_deg(i,:) = acosd(C2);

    [maxC1, id1] = max(C1);
    [maxC2, id2] = max(C2);

    K1 = K100 / maxC1;
    K2 = K110 / maxC2;

    K1_vals(i) = K1;
    K2_vals(i) = K2;

    if K1 <= K2
        minKeq_vals(i) = K1;
        best_family(i) = "{100}";
        best_plane_index(i) = id1;
    else
        minKeq_vals(i) = K2;
        best_family(i) = "{110}";
        best_plane_index(i) = id2;
    end

end

minKeq_MPa = minKeq_vals ./ 1e6;
K1_MPa = K1_vals ./ 1e6;
K2_MPa = K2_vals ./ 1e6;


figure('Color','w','Units','inches','Position',[0.5 0.5 5 4]); hold on;
histogram(minKeq_MPa, 40, ...
    'Normalization','pdf', ...
    'FaceColor',[0.65 0.65 0.65], ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$K_{\rm min}$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
ylabel('PDF','FontSize',16,'Interpreter','latex');
box on;
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex');

%% Finding K_T : Random Texture

clear; clc;

K100 = 1.78e6;
K110 = 1.78e6;

N = 1e6;

planes1 = [1 0 0; 0 1 0; 0 0 1]';
planes2 = [0 1 -1; 0 1 1; 1 -1 0; 1 0 -1; 1 0 1; 1 1 0]';
planes2 = planes2 ./ vecnorm(planes2);

C1_all = zeros(N,size(planes1,2));
C2_all = zeros(N,size(planes2,2));

K1_vals = zeros(N,1);
K2_vals = zeros(N,1);
minKeq_vals = zeros(N,1);

best_family = strings(N,1);
best_plane_index = zeros(N,1);

theta100_deg = zeros(N,size(planes1,2));
theta110_deg = zeros(N,size(planes2,2));

load_axis = [0 0 1];

for i = 1:N

    q = randn(4,1);
    q = q ./ norm(q);

    w  = q(1);
    x_ = q(2);
    y_ = q(3);
    z_ = q(4);

    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),    1-2*(x_^2+z_^2),   2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),    1-2*(x_^2+y_^2)];

    planes1_s = R * planes1;
    planes2_s = R * planes2;

    C1 = load_axis * planes1_s;
    C2 = load_axis * planes2_s;

    C1 = abs(C1);
    C2 = abs(C2);

    C1_all(i,:) = C1;
    C2_all(i,:) = C2;

    theta100_deg(i,:) = acosd(C1);
    theta110_deg(i,:) = acosd(C2);

    [maxC1, id1] = max(C1);
    [maxC2, id2] = max(C2);

    K1 = K100 / maxC1;
    K2 = K110 / maxC2;

    K1_vals(i) = K1;
    K2_vals(i) = K2;

    if K1 <= K2
        minKeq_vals(i) = K1;
        best_family(i) = "{100}";
        best_plane_index(i) = id1;
    else
        minKeq_vals(i) = K2;
        best_family(i) = "{110}";
        best_plane_index(i) = id2;
    end

end

minKeq_MPa = minKeq_vals ./ 1e6;
K1_MPa = K1_vals ./ 1e6;
K2_MPa = K2_vals ./ 1e6;

figure('Color','w','Units','inches','Position',[0.5 0.5 5 4]); hold on;
histogram(minKeq_MPa, 40, ...
    'Normalization','pdf', ...
    'FaceColor',[0.65 0.65 0.65], ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$K_{\rm min}$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
ylabel('PDF','FontSize',16,'Interpreter','latex');
box on;
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex');


%% Finding K_T: bootstarp EBSD

clear; clc;

K100 = 1.78e6;
K110 = 1.78e6;

N_boot = 1e5;
rng(1);

planes1 = [1 0 0; 0 1 0; 0 0 1]';
planes2 = [0 1 -1; 0 1 1; 1 -1 0; 1 0 -1; 1 0 1; 1 1 0]';
planes2 = planes2 ./ vecnorm(planes2);

X = load('EBSD_grain_orientation_X3d.mat');
X_3d_real = X.X_3d;

N_real = size(X_3d_real,1);
idxBoot = randi(N_real, N_boot, 1);
X_3d = X_3d_real(idxBoot,:);

N = size(X_3d,1);

C1_all = zeros(N,size(planes1,2));
C2_all = zeros(N,size(planes2,2));

K1_vals = zeros(N,1);
K2_vals = zeros(N,1);
minKeq_vals = zeros(N,1);

best_family = strings(N,1);
best_plane_index = zeros(N,1);

theta100_deg = zeros(N,size(planes1,2));
theta110_deg = zeros(N,size(planes2,2));

load_axis = [0 0 1];

for i = 1:N

    phi1 = X_3d(i,1);
    Phi  = X_3d(i,2);
    phi2 = X_3d(i,3);

    Rz1 = [cos(phi1) -sin(phi1) 0; sin(phi1) cos(phi1) 0; 0 0 1];
    Rx  = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
    Rz2 = [cos(phi2) -sin(phi2) 0; sin(phi2) cos(phi2) 0; 0 0 1];

    R = Rz1 * Rx * Rz2;

    planes1_s = R * planes1;
    planes2_s = R * planes2;

    C1 = abs(load_axis * planes1_s);
    C2 = abs(load_axis * planes2_s);

    C1_all(i,:) = C1;
    C2_all(i,:) = C2;

    theta100_deg(i,:) = acosd(C1);
    theta110_deg(i,:) = acosd(C2);

    [maxC1, id1] = max(C1);
    [maxC2, id2] = max(C2);

    K1 = K100 / maxC1;
    K2 = K110 / maxC2;

    K1_vals(i) = K1;
    K2_vals(i) = K2;

    if K1 <= K2
        minKeq_vals(i) = K1;
        best_family(i) = "{100}";
        best_plane_index(i) = id1;
    else
        minKeq_vals(i) = K2;
        best_family(i) = "{110}";
        best_plane_index(i) = id2;
    end

end

minKeq_MPa = minKeq_vals ./ 1e6;
K1_MPa = K1_vals ./ 1e6;
K2_MPa = K2_vals ./ 1e6;


figure('Color','w','Units','inches','Position',[0.5 0.5 5 4]); hold on;
histogram(minKeq_MPa, 60, ...
    'Normalization','pdf', ...
    'FaceColor',[0.65 0.65 0.65], ...
    'EdgeColor','k', ...
    'LineWidth',0.7, ...
    'FaceAlpha',0.8);
xlabel('$K_{\rm min}$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
ylabel('PDF','FontSize',16,'Interpreter','latex');
box on;
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex');

% fprintf('Real EBSD grains: %d\n', N_real);
% fprintf('Bootstrap samples: %d\n', N_boot);
% fprintf('Mean Kmin: %.4f MPa sqrt(m)\n', mean(minKeq_MPa,'omitnan'));
% fprintf('Std Kmin: %.4f MPa sqrt(m)\n', std(minKeq_MPa,'omitnan'));
% fprintf('CV Kmin: %.4f %%\n', 100*std(minKeq_MPa,'omitnan')/mean(minKeq_MPa,'omitnan'));


%%

clear; clc; close all;

rng(1);

N = 1e6;

K100 = 1.78e6;
K110 = 1.78e6;

KIC_GB_raw = readmatrix('AK_data.csv');
theta_IV_raw = readmatrix('angle_selected_rad.csv');
theta_G_raw = readmatrix('angle_all_rad.csv');

KIC_GB_raw = KIC_GB_raw(:);
theta_IV_raw = theta_IV_raw(:);
theta_G_raw = theta_G_raw(:);

KIC_GB_raw = KIC_GB_raw(isfinite(KIC_GB_raw));
theta_IV_raw = theta_IV_raw(isfinite(theta_IV_raw));
theta_G_raw = theta_G_raw(isfinite(theta_G_raw));

KIC_GB_raw = KIC_GB_raw(KIC_GB_raw > 0);
theta_IV_raw = theta_IV_raw(theta_IV_raw >= 0 & theta_IV_raw < pi/2);
theta_G_raw = theta_G_raw(theta_G_raw >= 0 & theta_G_raw < pi/2);

KIC_GB_samples_V = KIC_GB_raw(randi(numel(KIC_GB_raw), N, 1));
theta_IV_samples = theta_IV_raw(randi(numel(theta_IV_raw), N, 1));

K_V_samples = KIC_GB_samples_V ./ cos(theta_IV_samples);
K_V_samples = K_V_samples(isfinite(K_V_samples));
K_V_samples = K_V_samples(K_V_samples > 0);

KIC_GB_samples_G = KIC_GB_raw(randi(numel(KIC_GB_raw), N, 1));
theta_G_samples = theta_G_raw(randi(numel(theta_G_raw), N, 1));

K_G_samples = KIC_GB_samples_G ./ cos(theta_G_samples);
K_G_samples = K_G_samples(isfinite(K_G_samples));
K_G_samples = K_G_samples(K_G_samples > 0);

planes1 = [1 0 0; 0 1 0; 0 0 1]';
planes2 = [0 1 -1; 0 1 1; 1 -1 0; 1 0 -1; 1 0 1; 1 1 0]';
planes2 = planes2 ./ vecnorm(planes2);

X = load('EBSD_grain_orientation_X3d.mat');
X_3d_real = X.X_3d;

N_real = size(X_3d_real,1);
idxBoot = randi(N_real, N, 1);
X_3d = X_3d_real(idxBoot,:);

K_T_samples = zeros(N,1);
K100_samples = zeros(N,1);
K110_samples = zeros(N,1);

load_axis = [0 0 1];

for i = 1:N

    phi1 = X_3d(i,1);
    Phi  = X_3d(i,2);
    phi2 = X_3d(i,3);

    Rz1 = [cos(phi1) -sin(phi1) 0; sin(phi1) cos(phi1) 0; 0 0 1];
    Rx  = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
    Rz2 = [cos(phi2) -sin(phi2) 0; sin(phi2) cos(phi2) 0; 0 0 1];

    R = Rz1 * Rx * Rz2;

    planes1_s = R * planes1;
    planes2_s = R * planes2;

    C1 = abs(load_axis * planes1_s);
    C2 = abs(load_axis * planes2_s);

    maxC1 = max(C1);
    maxC2 = max(C2);

    K100_samples(i) = K100 / maxC1;
    K110_samples(i) = K110 / maxC2;
    K_T_samples(i) = min(K100_samples(i), K110_samples(i));

end

K_T_samples = K_T_samples ./ 1e6;
K100_samples = K100_samples ./ 1e6;
K110_samples = K110_samples ./ 1e6;

K_T_samples = K_T_samples(isfinite(K_T_samples));
K_T_samples = K_T_samples(K_T_samples > 0);

binEdges = 0:0.1:10;

K_V_plot = K_V_samples(K_V_samples >= binEdges(1) & K_V_samples <= binEdges(end));
K_G_plot = K_G_samples(K_G_samples >= binEdges(1) & K_G_samples <= binEdges(end));
K_T_plot = K_T_samples(K_T_samples >= binEdges(1) & K_T_samples <= binEdges(end));

figure('Color','w','Units','inches','Position',[0.5 0.5 12 3]);

subplot(1,3,1); hold on; box on; grid off;
histogram(K_V_plot, binEdges, ...
    'Normalization','pdf', ...
    'FaceColor',[0.1 0.35 0.85], ...
    'EdgeColor','k', ...
    'LineWidth',0.6, ...
    'FaceAlpha',0.75);
xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
ylabel('$f_K^{\rm I_V}$','FontSize',16,'Interpreter','latex');
xlim([1 5]);
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

subplot(1,3,2); hold on; box on; grid off;
histogram(K_G_plot, binEdges, ...
    'Normalization','pdf', ...
    'FaceColor','r', ...
    'EdgeColor','k', ...
    'LineWidth',0.6, ...
    'FaceAlpha',0.75);
xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
ylabel('$f_K^{\rm I}$','FontSize',16,'Interpreter','latex');
xlim([1 5]);
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

subplot(1,3,3); hold on; box on; grid off;
histogram(K_T_plot, 30, ...
    'Normalization','pdf', ...
    'FaceColor',[0.65 0.65 0.65], ...
    'EdgeColor','k', ...
    'LineWidth',0.6, ...
    'FaceAlpha',0.75);
xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
ylabel('$f_K^{\rm T}$','FontSize',16,'Interpreter','latex');
xlim([1.75 2.2]);
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');


Z1_samples = min(K_V_samples, K_T_samples);
Z2_samples = min(K_G_samples, K_T_samples); 

pV_empirical = mean(K_V_samples < K_T_samples);
pG_empirical = mean(K_G_samples < K_T_samples);

fprintf('Empirical P(K_V < K_T) = %.4f\n', pV_empirical);
fprintf('Empirical P(K_G < K_T) = %.4f\n', pG_empirical);
