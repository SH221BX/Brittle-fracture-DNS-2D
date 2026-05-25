clear; clc; close all;

rng(1);

N = 1e6;
n_segments = 100;
aVals = [1 1.125 1.25 1.5 1.75 2 2.25 2.5 2.75 3];

% KIC_GB_raw = readmatrix('AK_data.csv');                         % EBSD 
% theta_IV_raw = readmatrix('angle_selected_rad.csv');
% theta_G_raw = readmatrix('angle_all_rad.csv');


KIC_GB_raw = readmatrix('AK_data.csv');                          % Vornoi 3D
S_AT = load('A_T.mat');
S_AG = load('A_G.mat');

fnT = fieldnames(S_AT); fnG = fieldnames(S_AG);
theta_IV_raw = S_AT.(fnT{1}); theta_G_raw = S_AG.(fnG{1});
theta_IV_raw = deg2rad(theta_IV_raw(:));
theta_G_raw = deg2rad(theta_G_raw(:));


% KIC_GB_raw = readmatrix('AK_data.csv');                           % Voronoi 2D
% S_AT = load('T_Iv.mat'); S_AG = load('T_I.mat');

% fnT = fieldnames(S_AT); fnG = fieldnames(S_AG);
% theta_IV_raw = S_AT.(fnT{1}); theta_G_raw = S_AG.(fnG{1});
% theta_IV_raw = (theta_IV_raw(:)); theta_G_raw = (theta_G_raw(:));


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

% X = load('EBSD_grain_orientation_X3d.mat');
% X_3d_real = X.X_3d;

X = load('Random_texture_X3d.mat');
X_3d_real = X.X_3d_random;

N_real = size(X_3d_real,1);
Kgeom_real = zeros(N_real,1);
load_axis = [0 0 1];

for i = 1:N_real

    phi1 = X_3d_real(i,1);
    Phi  = X_3d_real(i,2);
    phi2 = X_3d_real(i,3);

    Rz1 = [cos(phi1) -sin(phi1) 0; sin(phi1) cos(phi1) 0; 0 0 1];
    Rx  = [1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
    Rz2 = [cos(phi2) -sin(phi2) 0; sin(phi2) cos(phi2) 0; 0 0 1];

    R = Rz1 * Rx * Rz2;

    C1 = abs(load_axis * (R * planes1));
    C2 = abs(load_axis * (R * planes2));

    Kgeom_real(i) = min(1 / max(C1), 1 / max(C2));

end

Kgeom_real = Kgeom_real(isfinite(Kgeom_real));
Kgeom_real = Kgeom_real(Kgeom_real > 0);

Nvec = 1:n_segments;
idxErr = [1 10:10:n_segments];
idxErr = idxErr(idxErr <= numel(Nvec));

nA = numel(aVals);

KM_mean_all = zeros(nA,n_segments);
KM_q25_all = zeros(nA,n_segments);
KM_q75_all = zeros(nA,n_segments);
KM_min_all = zeros(nA,n_segments);
KM_max_all = zeros(nA,n_segments);

KT_frac_mean_all = zeros(nA,n_segments);
KT_frac_q25_all = zeros(nA,n_segments);
KT_frac_q75_all = zeros(nA,n_segments);
KT_frac_min_all = zeros(nA,n_segments);
KT_frac_max_all = zeros(nA,n_segments);

P_T_all = zeros(nA,n_segments);
pV_all = zeros(nA,n_segments);

colors = [0.8500 0.1250 0.1250; 0.1216 0.4706 0.7059; 0.1000 0.5500 0.2000];

for aa = 1:nA

    a = aVals(aa);
    K_TX = a * mean(KIC_GB_raw,'omitnan');
    K_T_real = K_TX * Kgeom_real;
    K_T_real = K_T_real(isfinite(K_T_real));
    K_T_real = K_T_real(K_T_real > 0);
    K_T_samples = K_T_real(randi(numel(K_T_real), N, 1));

    pV = zeros(1,n_segments);
    P_T = zeros(1,n_segments);

    KM_mean = zeros(1,n_segments);
    KM_q25 = zeros(1,n_segments);
    KM_q75 = zeros(1,n_segments);
    KM_min = zeros(1,n_segments);
    KM_max = zeros(1,n_segments);

    KT_count_mean = zeros(1,n_segments);
    KT_count_q25 = zeros(1,n_segments);
    KT_count_q75 = zeros(1,n_segments);
    KT_count_min = zeros(1,n_segments);
    KT_count_max = zeros(1,n_segments);

    KV = K_V_samples(randi(numel(K_V_samples), N, 1));
    KT = K_T_samples(randi(numel(K_T_samples), N, 1));

    Z1fresh = min(KV, KT);
    isV = KV < KT;
    isT_first = ~isV;

    pV(1) = mean(isV);
    P_T(1) = 1 - pV(1);

    MaxSegment = Z1fresh;
    KT_count_path = double(isT_first);

    KM_mean(1) = mean(MaxSegment,'omitnan');
    KM_q25(1) = quantile(MaxSegment,0.25);
    KM_q75(1) = quantile(MaxSegment,0.75);
    KM_min(1) = min(MaxSegment);
    KM_max(1) = max(MaxSegment);

    KT_count_mean(1) = mean(KT_count_path);
    KT_count_q25(1) = quantile(KT_count_path,0.25);
    KT_count_q75(1) = quantile(KT_count_path,0.75);
    KT_count_min(1) = min(KT_count_path);
    KT_count_max(1) = max(KT_count_path);

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
        KT_count_path = KT_count_path + double(isT_Seg);

        KM_mean(seg) = mean(MaxSegment,'omitnan');
        KM_q25(seg) = quantile(MaxSegment,0.25);
        KM_q75(seg) = quantile(MaxSegment,0.75);
        KM_min(seg) = min(MaxSegment);
        KM_max(seg) = max(MaxSegment);

        KT_count_mean(seg) = mean(KT_count_path);
        KT_count_q25(seg) = quantile(KT_count_path,0.25);
        KT_count_q75(seg) = quantile(KT_count_path,0.75);
        KT_count_min(seg) = min(KT_count_path);
        KT_count_max(seg) = max(KT_count_path);

    end

    KM_mean_all(aa,:) = KM_mean;
    KM_q25_all(aa,:) = KM_q25;
    KM_q75_all(aa,:) = KM_q75;
    KM_min_all(aa,:) = KM_min;
    KM_max_all(aa,:) = KM_max;

    KT_frac_mean_all(aa,:) = KT_count_mean ./ Nvec;
    KT_frac_q25_all(aa,:) = KT_count_q25 ./ Nvec;
    KT_frac_q75_all(aa,:) = KT_count_q75 ./ Nvec;
    KT_frac_min_all(aa,:) = KT_count_min ./ Nvec;
    KT_frac_max_all(aa,:) = KT_count_max ./ Nvec;

    P_T_all(aa,:) = P_T;
    pV_all(aa,:) = pV;

    fprintf('a = %.1f\n', a);
    fprintf('K_TX = %.6f MPa sqrt(m)\n', K_TX);
    fprintf('Expected fraction of K_T segments = %.6f\n', KT_frac_mean_all(aa,end));
    fprintf('Final P_T = %.6f\n', P_T(end));
    fprintf('Final P_I = %.6f\n', pV(end));
    fprintf('E[K_M] final = %.6f MPa sqrt(m)\n\n', KM_mean(end));

end

% figure('Color','w','Units','inches','Position',[0.5 0.5 6 4]); hold on; box on; grid off;
% 
% for aa = 1:nA
%     plot(Nvec, KM_mean_all(aa,:), '-', 'Color',colors(aa,:), 'LineWidth',2.2, 'DisplayName',sprintf('$a = %.0f$',aVals(aa)));
%     plot(Nvec, KM_q25_all(aa,:), '--', 'Color',colors(aa,:), 'LineWidth',1.0, 'HandleVisibility','off');
%     plot(Nvec, KM_q75_all(aa,:), '--', 'Color',colors(aa,:), 'LineWidth',1.0, 'HandleVisibility','off');
% end
% 
% xlabel('$\rm{N}$','FontSize',16,'Interpreter','latex');
% ylabel('$\rm{K_R}$ (MPa$\sqrt{\rm m}$)','FontSize',16,'Interpreter','latex');
% legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% xlim([-5 105]);
% 
% figure('Color','w','Units','inches','Position',[0.5 0.5 6 4]); hold on; box on; grid off;
% 
% for aa = 1:nA
%     plot(Nvec, KT_frac_mean_all(aa,:), '-', 'Color',colors(aa,:), 'LineWidth',2.2, 'DisplayName',sprintf('$a = %.0f$',aVals(aa)));
%     plot(Nvec, KT_frac_q25_all(aa,:), '--', 'Color',colors(aa,:), 'LineWidth',1.0, 'HandleVisibility','off');
%     plot(Nvec, KT_frac_q75_all(aa,:), '--', 'Color',colors(aa,:), 'LineWidth',1.0, 'HandleVisibility','off');
% end
% 
% xlabel('$\rm{N}$','FontSize',16,'Interpreter','latex');
% ylabel('$\rm{N_T/N}$','FontSize',16,'Interpreter','latex');
% legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
% set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
% xlim([-5 105]);
% ylim([0 1]);

%%
a_plot = aVals(:);
NTN_plot = KT_frac_mean_all(:,end);
NTN_q25 = KT_frac_q25_all(:,end);
NTN_q75 = KT_frac_q75_all(:,end);
NTN_min = KT_frac_min_all(:,end);
NTN_max = KT_frac_max_all(:,end);

figure('Color','w','Units','inches','Position',[0.5 0.5 5 4]); hold on; box on; grid off;
errorbar(a_plot, NTN_plot, NTN_plot-NTN_min, NTN_max-NTN_plot, 'co', 'LineWidth',1.2, 'MarkerFaceColor','w', 'MarkerSize',6, 'CapSize',8, 'DisplayName','Range');
plot(a_plot, NTN_plot, 'c-', 'LineWidth',1.8, 'HandleVisibility','off');
plot(a_plot, NTN_q25, 'c--', 'LineWidth',1.0, 'DisplayName','IQR');
plot(a_plot, NTN_q75, 'c--', 'LineWidth',1.0, 'HandleVisibility','off');
scatter(a_plot, NTN_plot, 45, 'k', 'filled', 'HandleVisibility','off');
xlabel('$\rm{R}$','FontSize',16,'Interpreter','latex');
ylabel('$\mathrm{\kappa}$','FontSize',16,'Interpreter','latex');
legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
xlim([min(a_plot)-0.2 max(a_plot)+0.2]);
ylim([0 1]);


%%

a_plot = aVals(:);
NTN_plot = KT_frac_mean_all(:,end);
NTN_q25 = KT_frac_q25_all(:,end);
NTN_q75 = KT_frac_q75_all(:,end);
NTN_min = KT_frac_min_all(:,end);
NTN_max = KT_frac_max_all(:,end);

NTN_plot(NTN_plot<=0) = NaN;
NTN_q25(NTN_q25<=0) = NaN;
NTN_q75(NTN_q75<=0) = NaN;
NTN_min(NTN_min<=0) = NaN;
NTN_max(NTN_max<=0) = NaN;

figure('Color','w','Units','inches','Position',[0.5 0.5 5 4]); hold on; box on; grid off;
% errorbar(a_plot, NTN_plot, NTN_plot-NTN_min, NTN_max-NTN_plot, 'ko', 'LineWidth',1.2, 'MarkerFaceColor','w', 'MarkerSize',6, 'CapSize',8, 'DisplayName','Range');
plot(a_plot, NTN_plot, 'k-', 'LineWidth',1.8, 'HandleVisibility','off');
% plot(a_plot, NTN_q25, 'k--', 'LineWidth',1.0, 'DisplayName','IQR');
% plot(a_plot, NTN_q75, 'k--', 'LineWidth',1.0, 'HandleVisibility','off');
scatter(a_plot, NTN_plot, 45, 'k', 'filled', 'HandleVisibility','off');
xlabel('$\rm{R}$','FontSize',16,'Interpreter','latex');
ylabel('$\left.\mathrm{N_T/N}\right|_{N=100}$','FontSize',16,'Interpreter','latex');
legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman','XScale','log','YScale','log');
xlim([min(a_plot)*0.9 max(a_plot)*1.1]);


