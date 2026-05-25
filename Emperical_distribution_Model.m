%% EBSD Texture/ MEAM GB Library

clear; clc; 

rng(1);
N = 1e6; n_segments = 100;
K100 = 1.78e6; K110 = 1.78e6;

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

K_T_real = zeros(N_real,1);
K100_real = zeros(N_real,1);
K110_real = zeros(N_real,1);

load_axis = [0 0 1];

for i = 1:N_real

    phi1 = X_3d_real(i,1);
    Phi  = X_3d_real(i,2);
    phi2 = X_3d_real(i,3);

    Rz1 = [cos(phi1) -sin(phi1) 0; sin(phi1) cos(phi1) 0; 0 0 1];
    Rx  = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
    Rz2 = [cos(phi2) -sin(phi2) 0; sin(phi2) cos(phi2) 0; 0 0 1];

    R = Rz1 * Rx * Rz2;

    planes1_s = R * planes1;
    planes2_s = R * planes2;

    C1 = abs(load_axis * planes1_s);
    C2 = abs(load_axis * planes2_s);

    K100_real(i) = K100 / max(C1);
    K110_real(i) = K110 / max(C2);
    K_T_real(i) = min(K100_real(i), K110_real(i));

end

K_T_real = K_T_real ./ 1e6;
K100_real = K100_real ./ 1e6;
K110_real = K110_real ./ 1e6;

K_T_real = K_T_real(isfinite(K_T_real));
K_T_real = K_T_real(K_T_real > 0);

K_T_samples = K_T_real(randi(numel(K_T_real), N, 1));

binEdges = 0:0.1:10;

K_V_plot = K_V_samples(K_V_samples >= binEdges(1) & K_V_samples <= binEdges(end));
K_G_plot = K_G_samples(K_G_samples >= binEdges(1) & K_G_samples <= binEdges(end));
K_T_plot = K_T_samples(K_T_samples >= binEdges(1) & K_T_samples <= binEdges(end));

% figure('Color','w','Units','inches','Position',[0.5 0.5 12 3]);
% 
% subplot(1,3,1); hold on; box on; grid off;
% histogram(K_V_plot, binEdges, ...
%     'Normalization','pdf', ...
%     'FaceColor',[0.1 0.35 0.85], ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm I_V}$','FontSize',16,'Interpreter','latex');
% xlim([1 5]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
% subplot(1,3,2); hold on; box on; grid off;
% histogram(K_G_plot, binEdges, ...
%     'Normalization','pdf', ...
%     'FaceColor','r', ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm I}$','FontSize',16,'Interpreter','latex');
% xlim([1 5]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
% subplot(1,3,3); hold on; box on; grid off;
% histogram(K_T_plot, 30, ...
%     'Normalization','pdf', ...
%     'FaceColor',[0.65 0.65 0.65], ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm T}$','FontSize',16,'Interpreter','latex');
% xlim([1.75 2.2]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

pV = zeros(1,n_segments);
P_T = zeros(1,n_segments);
mean_KM = zeros(1,n_segments);

KV = K_V_samples(randi(numel(K_V_samples), N, 1));
KT = K_T_samples(randi(numel(K_T_samples), N, 1));

Z1fresh = min(KV, KT);
isV = KV < KT;

pV(1) = mean(isV);
P_T(1) = 1 - pV(1);

MaxSegment = Z1fresh;
mean_KM(1) = mean(MaxSegment,'omitnan');
KM_mean = zeros(1,n_segments);
KM_q25 = zeros(1,n_segments);
KM_q75 = zeros(1,n_segments);
KM_min = zeros(1,n_segments);
KM_max = zeros(1,n_segments);

KM_mean(1) = mean(MaxSegment,'omitnan');
KM_q25(1) = quantile(MaxSegment,0.25);
KM_q75(1) = quantile(MaxSegment,0.75);
KM_min(1) = min(MaxSegment);
KM_max(1) = max(MaxSegment);

for seg = 2:n_segments

    KV = K_V_samples(randi(numel(K_V_samples), N, 1));
    KG = K_G_samples(randi(numel(K_G_samples), N, 1));
    KT1 = K_T_samples(randi(numel(K_T_samples), N, 1));
    KT2 = K_T_samples(randi(numel(K_T_samples), N, 1));

    Z1fresh = min(KV, KT1);
    Z2fresh = min(KG, KT2);

    isV = KV < KT1;
    isG = KG < KT2;

    mixFlag = rand(N,1) <= pV(seg-1);

    Segment = zeros(N,1);
    Segment(mixFlag) = Z1fresh(mixFlag);
    Segment(~mixFlag) = Z2fresh(~mixFlag);

    isT_Seg = false(N,1);
    isT_Seg(mixFlag) = ~isV(mixFlag);
    isT_Seg(~mixFlag) = ~isG(~mixFlag);

    P_T(seg) = mean(isT_Seg);
    pV(seg) = 1 - P_T(seg);

    MaxSegment = max(MaxSegment, Segment);

    mean_KM(seg) = mean(MaxSegment,'omitnan');

    KM_mean(seg) = mean(MaxSegment,'omitnan');
    KM_q25(seg) = quantile(MaxSegment,0.25);
    KM_q75(seg) = quantile(MaxSegment,0.75);
    KM_min(seg) = min(MaxSegment);
    KM_max(seg) = max(MaxSegment);

end

Nvec = 1:n_segments;
idxErr = [1 10:10:100];
idxErr = idxErr(idxErr <= numel(Nvec));

figure('Color','w','Units','inches','Position',[0.5 0.5 6 4]); hold on; box on; grid off;
fill([Nvec fliplr(Nvec)], [KM_q25 fliplr(KM_q75)], 'r', 'EdgeColor','none', 'FaceAlpha',0.25, 'HandleVisibility','off');
plot(Nvec, KM_mean, 'r-', 'LineWidth',2.2, 'DisplayName','$\rm{E[K_M]}$');
plot(Nvec, KM_q25, 'r--', 'LineWidth',1.2, 'DisplayName','IQR');
plot(Nvec, KM_q75, 'r--', 'LineWidth',1.2, 'HandleVisibility','off');

errorbar(Nvec(idxErr), KM_mean(idxErr), KM_mean(idxErr)-KM_min(idxErr), KM_max(idxErr)-KM_mean(idxErr), ...
    'ro', 'LineWidth',1.2, 'MarkerFaceColor','w', 'MarkerSize',5, 'CapSize',7, 'DisplayName','Range');

xlabel('$\rm{N}$','FontSize',16,'Interpreter','latex');
ylabel('$\rm{K_R}$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
xlim([-5 105]);

fprintf('P(K_V < K_T) initial = %.4f\n', pV(1));
fprintf('P_T initial = %.4f\n', P_T(1));
fprintf('P_I final = %.4f\n', pV(end));
fprintf('P_T final = %.4f\n', P_T(end));
fprintf('E[K_M] final = %.4f MPa sqrt(m)\n', mean_KM(end));

% figure('Color','w','Units','inches','Position',[0.5 0.5 6 4]); hold on; box on; grid off;
% fill([Nvec fliplr(Nvec)], [KM_min fliplr(KM_max)], [0.75 0.75 0.75], 'EdgeColor','none', 'FaceAlpha',0.25, 'DisplayName','Min--max');
% fill([Nvec fliplr(Nvec)], [KM_q25 fliplr(KM_q75)], [0.1216 0.4706 0.7059], 'EdgeColor','none', 'FaceAlpha',0.25, 'DisplayName','IQR');
% plot(Nvec, KM_mean, 'b-', 'LineWidth',2.2, 'DisplayName','$E[K_M]$');
% plot(Nvec, KM_q25, 'b--', 'LineWidth',1.2, 'HandleVisibility','off');
% plot(Nvec, KM_q75, 'b--', 'LineWidth',1.2, 'HandleVisibility','off');

% figure('Color','w','Units','inches','Position',[0.5 0.5 5 4]); hold on; box on; grid off;
% plot(1:n_segments, pV, 'r-', 'LineWidth',2.2); hold on;
% plot(1:n_segments, P_T, 'k-', 'LineWidth',2.2);
% xlabel('$N$','FontSize',16,'Interpreter','latex');
% ylabel('Probability','FontSize',16,'Interpreter','latex');
% legend({'$P_I$','$P_T$'},'Interpreter','latex','Location','best','FontSize',13);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');









%% Voronoi 3D microstructure/ MEAM GB Library

clear; clc; 

rng(1);
N = 1e6; n_segments = 100;
K100 = 1.78e6; K110 = 1.78e6;

KIC_GB_raw = readmatrix('AK_data.csv');

S_AT = load('A_T.mat');
S_AG = load('A_G.mat');

fnT = fieldnames(S_AT);
fnG = fieldnames(S_AG);

theta_IV_raw = S_AT.(fnT{1});
theta_G_raw = S_AG.(fnG{1});

theta_IV_raw = deg2rad(theta_IV_raw(:));
theta_G_raw = deg2rad(theta_G_raw(:));


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

K_T_real = zeros(N_real,1);
K100_real = zeros(N_real,1);
K110_real = zeros(N_real,1);

load_axis = [0 0 1];

for i = 1:N_real

    phi1 = X_3d_real(i,1);
    Phi  = X_3d_real(i,2);
    phi2 = X_3d_real(i,3);

    Rz1 = [cos(phi1) -sin(phi1) 0; sin(phi1) cos(phi1) 0; 0 0 1];
    Rx  = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
    Rz2 = [cos(phi2) -sin(phi2) 0; sin(phi2) cos(phi2) 0; 0 0 1];

    R = Rz1 * Rx * Rz2;

    planes1_s = R * planes1;
    planes2_s = R * planes2;

    C1 = abs(load_axis * planes1_s);
    C2 = abs(load_axis * planes2_s);

    K100_real(i) = K100 / max(C1);
    K110_real(i) = K110 / max(C2);
    K_T_real(i) = min(K100_real(i), K110_real(i));

end

K_T_real = K_T_real ./ 1e6;
K100_real = K100_real ./ 1e6;
K110_real = K110_real ./ 1e6;

K_T_real = K_T_real(isfinite(K_T_real));
K_T_real = K_T_real(K_T_real > 0);

K_T_samples = K_T_real(randi(numel(K_T_real), N, 1));

binEdges = 0:0.1:10;

K_V_plot = K_V_samples(K_V_samples >= binEdges(1) & K_V_samples <= binEdges(end));
K_G_plot = K_G_samples(K_G_samples >= binEdges(1) & K_G_samples <= binEdges(end));
K_T_plot = K_T_samples(K_T_samples >= binEdges(1) & K_T_samples <= binEdges(end));

% figure('Color','w','Units','inches','Position',[0.5 0.5 12 3]);
% 
% subplot(1,3,1); hold on; box on; grid off;
% histogram(K_V_plot, binEdges, ...
%     'Normalization','pdf', ...
%     'FaceColor',[0.1 0.35 0.85], ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm I_V}$','FontSize',16,'Interpreter','latex');
% xlim([1 5]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
% subplot(1,3,2); hold on; box on; grid off;
% histogram(K_G_plot, binEdges, ...
%     'Normalization','pdf', ...
%     'FaceColor','r', ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm I}$','FontSize',16,'Interpreter','latex');
% xlim([1 5]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
% subplot(1,3,3); hold on; box on; grid off;
% histogram(K_T_plot, 30, ...
%     'Normalization','pdf', ...
%     'FaceColor',[0.65 0.65 0.65], ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm T}$','FontSize',16,'Interpreter','latex');
% xlim([1.75 2.2]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

pV = zeros(1,n_segments);
P_T = zeros(1,n_segments);
mean_KM = zeros(1,n_segments);

KV = K_V_samples(randi(numel(K_V_samples), N, 1));
KT = K_T_samples(randi(numel(K_T_samples), N, 1));

Z1fresh = min(KV, KT);
isV = KV < KT;

pV(1) = mean(isV);
P_T(1) = 1 - pV(1);

MaxSegment = Z1fresh;
mean_KM(1) = mean(MaxSegment,'omitnan');
KM_mean = zeros(1,n_segments);
KM_q25 = zeros(1,n_segments);
KM_q75 = zeros(1,n_segments);
KM_min = zeros(1,n_segments);
KM_max = zeros(1,n_segments);

KM_mean(1) = mean(MaxSegment,'omitnan');
KM_q25(1) = quantile(MaxSegment,0.25);
KM_q75(1) = quantile(MaxSegment,0.75);
KM_min(1) = min(MaxSegment);
KM_max(1) = max(MaxSegment);

for seg = 2:n_segments

    KV = K_V_samples(randi(numel(K_V_samples), N, 1));
    KG = K_G_samples(randi(numel(K_G_samples), N, 1));
    KT1 = K_T_samples(randi(numel(K_T_samples), N, 1));
    KT2 = K_T_samples(randi(numel(K_T_samples), N, 1));

    Z1fresh = min(KV, KT1);
    Z2fresh = min(KG, KT2);

    isV = KV < KT1;
    isG = KG < KT2;

    mixFlag = rand(N,1) <= pV(seg-1);

    Segment = zeros(N,1);
    Segment(mixFlag) = Z1fresh(mixFlag);
    Segment(~mixFlag) = Z2fresh(~mixFlag);

    isT_Seg = false(N,1);
    isT_Seg(mixFlag) = ~isV(mixFlag);
    isT_Seg(~mixFlag) = ~isG(~mixFlag);

    P_T(seg) = mean(isT_Seg);
    pV(seg) = 1 - P_T(seg);

    MaxSegment = max(MaxSegment, Segment);

    mean_KM(seg) = mean(MaxSegment,'omitnan');

    KM_mean(seg) = mean(MaxSegment,'omitnan');
    KM_q25(seg) = quantile(MaxSegment,0.25);
    KM_q75(seg) = quantile(MaxSegment,0.75);
    KM_min(seg) = min(MaxSegment);
    KM_max(seg) = max(MaxSegment);

end

Nvec = 1:n_segments;
idxErr = [1 10:10:100];
idxErr = idxErr(idxErr <= numel(Nvec));

figure('Color','w','Units','inches','Position',[0.5 0.5 6 4]); hold on; box on; grid off;
fill([Nvec fliplr(Nvec)], [KM_q25 fliplr(KM_q75)], 'b', 'EdgeColor','none', 'FaceAlpha',0.25, 'HandleVisibility','off');
plot(Nvec, KM_mean, 'b-', 'LineWidth',2.2, 'DisplayName','$\rm{E[K_M]}$');
plot(Nvec, KM_q25, 'b--', 'LineWidth',1.2, 'DisplayName','IQR');
plot(Nvec, KM_q75, 'b--', 'LineWidth',1.2, 'HandleVisibility','off');

errorbar(Nvec(idxErr), KM_mean(idxErr), KM_mean(idxErr)-KM_min(idxErr), KM_max(idxErr)-KM_mean(idxErr), ...
    'bo', 'LineWidth',1.2, 'MarkerFaceColor','w', 'MarkerSize',5, 'CapSize',7, 'DisplayName','Range');

xlabel('$\rm{N}$','FontSize',16,'Interpreter','latex');
ylabel('$\rm{K_R}$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
xlim([-5 105]);

fprintf('P(K_V < K_T) initial = %.4f\n', pV(1));
fprintf('P_T initial = %.4f\n', P_T(1));
fprintf('P_I final = %.4f\n', pV(end));
fprintf('P_T final = %.4f\n', P_T(end));
fprintf('E[K_M] final = %.4f MPa sqrt(m)\n', mean_KM(end));





%% Voronoi 2D microstructure/ MEAM GB Library

clear; clc; 

rng(1);
N = 1e6; n_segments = 100;
K100 = 1.78e6; K110 = 1.78e6;

KIC_GB_raw = readmatrix('AK_data.csv');
S_AT = load('T_Iv.mat'); S_AG = load('T_I.mat');

fnT = fieldnames(S_AT); fnG = fieldnames(S_AG);
theta_IV_raw = S_AT.(fnT{1}); theta_G_raw = S_AG.(fnG{1});
theta_IV_raw = (theta_IV_raw(:)); theta_G_raw = (theta_G_raw(:));


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

K_T_real = zeros(N_real,1);
K100_real = zeros(N_real,1);
K110_real = zeros(N_real,1);

load_axis = [0 0 1];

for i = 1:N_real

    phi1 = X_3d_real(i,1);
    Phi  = X_3d_real(i,2);
    phi2 = X_3d_real(i,3);

    Rz1 = [cos(phi1) -sin(phi1) 0; sin(phi1) cos(phi1) 0; 0 0 1];
    Rx  = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
    Rz2 = [cos(phi2) -sin(phi2) 0; sin(phi2) cos(phi2) 0; 0 0 1];

    R = Rz1 * Rx * Rz2;

    planes1_s = R * planes1;
    planes2_s = R * planes2;

    C1 = abs(load_axis * planes1_s);
    C2 = abs(load_axis * planes2_s);

    K100_real(i) = K100 / max(C1);
    K110_real(i) = K110 / max(C2);
    K_T_real(i) = min(K100_real(i), K110_real(i));

end

K_T_real = K_T_real ./ 1e6;
K100_real = K100_real ./ 1e6;
K110_real = K110_real ./ 1e6;

K_T_real = K_T_real(isfinite(K_T_real));
K_T_real = K_T_real(K_T_real > 0);

K_T_samples = K_T_real(randi(numel(K_T_real), N, 1));

binEdges = 0:0.1:10;

K_V_plot = K_V_samples(K_V_samples >= binEdges(1) & K_V_samples <= binEdges(end));
K_G_plot = K_G_samples(K_G_samples >= binEdges(1) & K_G_samples <= binEdges(end));
K_T_plot = K_T_samples(K_T_samples >= binEdges(1) & K_T_samples <= binEdges(end));

% figure('Color','w','Units','inches','Position',[0.5 0.5 12 3]);
% 
% subplot(1,3,1); hold on; box on; grid off;
% histogram(K_V_plot, binEdges, ...
%     'Normalization','pdf', ...
%     'FaceColor',[0.1 0.35 0.85], ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm I_V}$','FontSize',16,'Interpreter','latex');
% xlim([1 5]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
% subplot(1,3,2); hold on; box on; grid off;
% histogram(K_G_plot, binEdges, ...
%     'Normalization','pdf', ...
%     'FaceColor','r', ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm I}$','FontSize',16,'Interpreter','latex');
% xlim([1 5]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
% subplot(1,3,3); hold on; box on; grid off;
% histogram(K_T_plot, 30, ...
%     'Normalization','pdf', ...
%     'FaceColor',[0.65 0.65 0.65], ...
%     'EdgeColor','k', ...
%     'LineWidth',0.6, ...
%     'FaceAlpha',0.75);
% xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% ylabel('$f_K^{\rm T}$','FontSize',16,'Interpreter','latex');
% xlim([1.75 2.2]);
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');

pV = zeros(1,n_segments);
P_T = zeros(1,n_segments);
mean_KM = zeros(1,n_segments);

KV = K_V_samples(randi(numel(K_V_samples), N, 1));
KT = K_T_samples(randi(numel(K_T_samples), N, 1));

Z1fresh = min(KV, KT);
isV = KV < KT;

pV(1) = mean(isV);
P_T(1) = 1 - pV(1);

MaxSegment = Z1fresh;
mean_KM(1) = mean(MaxSegment,'omitnan');
KM_mean = zeros(1,n_segments);
KM_q25 = zeros(1,n_segments);
KM_q75 = zeros(1,n_segments);
KM_min = zeros(1,n_segments);
KM_max = zeros(1,n_segments);

KM_mean(1) = mean(MaxSegment,'omitnan');
KM_q25(1) = quantile(MaxSegment,0.25);
KM_q75(1) = quantile(MaxSegment,0.75);
KM_min(1) = min(MaxSegment);
KM_max(1) = max(MaxSegment);

for seg = 2:n_segments

    KV = K_V_samples(randi(numel(K_V_samples), N, 1));
    KG = K_G_samples(randi(numel(K_G_samples), N, 1));
    KT1 = K_T_samples(randi(numel(K_T_samples), N, 1));
    KT2 = K_T_samples(randi(numel(K_T_samples), N, 1));

    Z1fresh = min(KV, KT1);
    Z2fresh = min(KG, KT2);

    isV = KV < KT1;
    isG = KG < KT2;

    mixFlag = rand(N,1) <= pV(seg-1);

    Segment = zeros(N,1);
    Segment(mixFlag) = Z1fresh(mixFlag);
    Segment(~mixFlag) = Z2fresh(~mixFlag);

    isT_Seg = false(N,1);
    isT_Seg(mixFlag) = ~isV(mixFlag);
    isT_Seg(~mixFlag) = ~isG(~mixFlag);

    P_T(seg) = mean(isT_Seg);
    pV(seg) = 1 - P_T(seg);

    MaxSegment = max(MaxSegment, Segment);

    mean_KM(seg) = mean(MaxSegment,'omitnan');

    KM_mean(seg) = mean(MaxSegment,'omitnan');
    KM_q25(seg) = quantile(MaxSegment,0.25);
    KM_q75(seg) = quantile(MaxSegment,0.75);
    KM_min(seg) = min(MaxSegment);
    KM_max(seg) = max(MaxSegment);

end

Nvec = 1:n_segments;
idxErr = [1 10:10:100];
idxErr = idxErr(idxErr <= numel(Nvec));

figure('Color','w','Units','inches','Position',[0.5 0.5 6 4]); hold on; box on; grid off;
fill([Nvec fliplr(Nvec)], [KM_q25 fliplr(KM_q75)], 'k', 'EdgeColor','none', 'FaceAlpha',0.25, 'HandleVisibility','off');
plot(Nvec, KM_mean, 'k-', 'LineWidth',2.2, 'DisplayName','$\rm{E[K_M]}$');
plot(Nvec, KM_q25, 'k--', 'LineWidth',1.2, 'DisplayName','IQR');
plot(Nvec, KM_q75, 'k--', 'LineWidth',1.2, 'HandleVisibility','off');

errorbar(Nvec(idxErr), KM_mean(idxErr), KM_mean(idxErr)-KM_min(idxErr), KM_max(idxErr)-KM_mean(idxErr), ...
    'ko', 'LineWidth',1.2, 'MarkerFaceColor','w', 'MarkerSize',5, 'CapSize',7, 'DisplayName','Range');

xlabel('$\rm{N}$','FontSize',16,'Interpreter','latex');
ylabel('$\rm{K_R}$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
xlim([-5 105]);

fprintf('P(K_V < K_T) initial = %.4f\n', pV(1));
fprintf('P_T initial = %.4f\n', P_T(1));
fprintf('P_I final = %.4f\n', pV(end));
fprintf('P_T final = %.4f\n', P_T(end));
fprintf('E[K_M] final = %.4f MPa sqrt(m)\n', mean_KM(end));









%%  Kmax evolving with N visualization
% 
% pV = zeros(1,n_segments);
% P_T = zeros(1,n_segments);
% mean_KM = zeros(1,n_segments);
% 
% binEdges_KM = linspace(1, 2.2, 100);
% 
% figHist = figure('Color','w','Units','inches','Position',[0.5 0.5 10 4]);
% 
% KV = K_V_samples(randi(numel(K_V_samples), N, 1));
% KT = K_T_samples(randi(numel(K_T_samples), N, 1));
% 
% Z1fresh = min(KV, KT);
% isV = KV < KT;
% 
% pV(1) = mean(isV);
% P_T(1) = 1 - pV(1);
% 
% Segment = Z1fresh;
% MaxSegment = Segment;
% mean_KM(1) = mean(MaxSegment,'omitnan');
% 
% for seg = 1:n_segments
% 
%     if seg > 1
% 
%         KV = K_V_samples(randi(numel(K_V_samples), N, 1));
%         KG = K_G_samples(randi(numel(K_G_samples), N, 1));
%         KT1 = K_T_samples(randi(numel(K_T_samples), N, 1));
%         KT2 = K_T_samples(randi(numel(K_T_samples), N, 1));
% 
%         Z1fresh = min(KV, KT1);
%         Z2fresh = min(KG, KT2);
% 
%         isV = KV < KT1;
%         isG = KG < KT2;
% 
%         mixFlag = rand(N,1) <= pV(seg-1);
% 
%         Segment = zeros(N,1);
%         Segment(mixFlag) = Z1fresh(mixFlag);
%         Segment(~mixFlag) = Z2fresh(~mixFlag);
% 
%         isT_Seg = false(N,1);
%         isT_Seg(mixFlag) = ~isV(mixFlag);
%         isT_Seg(~mixFlag) = ~isG(~mixFlag);
% 
%         P_T(seg) = mean(isT_Seg);
%         pV(seg) = 1 - P_T(seg);
% 
%         MaxSegment = max(MaxSegment, Segment);
%         mean_KM(seg) = mean(MaxSegment,'omitnan');
% 
%     end
% 
%     figure(figHist);
%     clf;
% 
%     % subplot(1,2,1); hold on; box on; grid off;
%     % histogram(Segment, binEdges_KM, ...
%     %     'Normalization','pdf', ...
%     %     'FaceColor',[0.85 0.15 0.15], ...
%     %     'EdgeColor','k', ...
%     %     'LineWidth',0.5, ...
%     %     'FaceAlpha',0.75);
%     % xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
%     % ylabel('$f_{Z_N}$','FontSize',16,'Interpreter','latex');
%     % title(['$Z_N,\ N = ' num2str(seg) '$'], ...
%     %     'FontSize',16, ...
%     %     'FontWeight','normal', ...
%     %     'Interpreter','latex');
%     % xlim([1 2.2]);
%     % set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
%     % subplot(1,2,2); hold on; box on; grid off;
%     histogram(MaxSegment, binEdges_KM,'Normalization','pdf', ...
%         'FaceColor',[0.7 0.4 0.0], ...
%         'EdgeColor','k', ...
%         'LineWidth',0.5, ...
%         'FaceAlpha',0.75);
%     xline(mean_KM(seg), 'k-', 'LineWidth',2.0);
%     xlabel('$K$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
%     ylabel('$f_{K_M}$','FontSize',16,'Interpreter','latex');
%     title(['$K_M(N),\ N = ' num2str(seg) '$'], ...
%         'FontSize',16, ...
%         'FontWeight','normal', ...
%         'Interpreter','latex');
%     xlim([1 2.2]);
%     set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% 
%     pause(0.5);
%     drawnow;
% 
% end
% 

