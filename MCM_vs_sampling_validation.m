clear all; clc; close all; warning off;

rng(1);

load('AK_data.mat');

K_GB_mean = mean(AK) * 1e6;
R0_lo = min(AK) * 1e6;
K_T_vals = linspace(K_GB_mean, 3*K_GB_mean, 10);

N_R0 = 2000;
X = 10000;
n_segments = 100;
a_beta = 1.9896;
b_beta = 1.5718;
r = 1.0;
phi0 = atan(1/sqrt(2));

R0_dist = AK(randi(numel(AK), N_R0, 1)) * 1e6;

theta1 = linspace(0, atan(2*sqrt(r)-sqrt(2)), X);
f1 = zeros(size(theta1));
t11 = atan(sqrt(2*r)-1);
i = theta1 <= t11;
f1(i) = sin(theta1(i));
c = (sqrt(2*r)-1).*cot(theta1(~i));
f1(~i) = (2/pi)*(asin(c)-acos(c)).*sin(theta1(~i));
n1 = trapz(theta1, f1);

u1 = linspace(1.0000001, sqrt(3-4*sqrt(2*r)+4*r), X);
th1 = acos(1./u1);
g1 = sin(th1);
i = th1 > t11;
c = (sqrt(2*r)-1).*cot(th1(i));
g1(i) = (2/pi)*(asin(c)-acos(c)).*sin(th1(i));
g1 = g1./(u1.*sqrt(u1.^2-1));
m1 = trapz(u1,g1);

t21 = atan(sqrt(2/r)-1);
t23 = atan(sqrt(3-4*sqrt(2*r)+3*r)/sqrt(r));
theta2 = linspace(0,t23,X);
f2 = zeros(size(theta2));
i1 = theta2 <= t21;
f2(i1) = sin(theta2(i1));
i2 = theta2 > t21 & theta2 <= pi/6;
c = (sqrt(2/r)-1).*cot(theta2(i2));
f2(i2) = (2/(pi-2*phi0))*(asin(c)-phi0).*sin(theta2(i2));
i3 = theta2 > pi/6 & theta2 <= t23;
c1 = (sqrt(2/r)-1).*cot(theta2(i3));
c2 = cot(theta2(i3))/sqrt(3);
f2(i3) = (2/(pi-2*phi0))*(asin(c1)-acos(c2)-phi0).*sin(theta2(i3));
n2 = trapz(theta2,f2);

u2 = linspace(sqrt(r),sqrt(3-4*sqrt(2*r)+4*r),X);
th2 = acos(sqrt(r)./u2);
g2 = zeros(size(th2));
j1 = th2 <= t21;
g2(j1) = sin(th2(j1));
j2 = th2 > t21 & th2 <= pi/6;
c = (sqrt(2/r)-1).*cot(th2(j2));
g2(j2) = (2/(pi-2*phi0))*(asin(c)-phi0).*sin(th2(j2));
j3 = th2 > pi/6 & th2 <= t23;
c1 = (sqrt(2/r)-1).*cot(th2(j3));
c2 = cot(th2(j3))/sqrt(3);
g2(j3) = (2/(pi-2*phi0))*(asin(c1)-acos(c2)-phi0).*sin(th2(j3));
g2 = g2.*(sqrt(r)./(u2.*sqrt(u2.^2-r)));
m2 = n2;

theta3 = linspace(0, phi0, X);
f3 = zeros(size(theta3));
i = theta3 <= pi/6;
f3(i) = sin(theta3(i));
c2 = cot(theta3(~i))/sqrt(3);
f3(~i) = (1/phi0)*(phi0-acos(c2)).*sin(theta3(~i));
n3 = trapz(theta3, f3);

u3 = linspace(sqrt(r), sqrt(r)/cos(phi0), X);
th3 = acos(sqrt(r)./u3);
g3 = sin(th3);
i = th3 > pi/6;
c2 = cot(th3(i))/sqrt(3);
g3(i) = (1/phi0)*(phi0-acos(c2)).*sin(th3(i));
g3 = g3.*(sqrt(r)./(u3.*sqrt(u3.^2-r)));
m3 = n3;

K_norm = linspace(1, sqrt(r)/cos(phi0), X);
rho_norm = ((pi/2)*interp1(u1,g1,K_norm,'linear',0) + (pi-2*phi0)*interp1(u2,g2,K_norm,'linear',0) + 2*phi0*interp1(u3,g3,K_norm,'linear',0)) / ((pi/2)*m1 + (pi-2*phi0)*m2 + 2*phi0*m3);

w_values = zeros(1,numel(K_T_vals));
y_terminal = zeros(1,numel(K_T_vals));
y_cumulative = zeros(1,numel(K_T_vals));
mean_pdf_all = zeros(numel(K_T_vals),n_segments);
pI_from_V_all = zeros(1,numel(K_T_vals));
pT_from_V_all = zeros(1,numel(K_T_vals));
pI_from_G_all = zeros(1,numel(K_T_vals));
pT_from_G_all = zeros(1,numel(K_T_vals));

for idx = 1:numel(K_T_vals)

    K_T = K_T_vals(idx);
    w = K_T / K_GB_mean;

    K = K_norm * K_T;
    rho = rho_norm / K_T;

    z_max = max(K);
    z = linspace(R0_lo + 1.0, z_max, X);

    fKV_mix_raw = zeros(1,numel(z));
    fKG_mix_raw = zeros(1,numel(z));

    for j = 1:N_R0

        R0_j = R0_dist(j);

        if R0_j >= z_max
            continue
        end

        fKV_j = @(x) 1./((pi/2)*beta(a_beta,b_beta)) .* ((acos(R0_j./x)/(pi/2)).^(a_beta-1)) .* (1 - acos(R0_j./x)/(pi/2)).^(b_beta-1) .* abs(R0_j./(x.^2.*sqrt(1-(R0_j./x).^2))) .* (x>=R0_j & x<=z_max);
        fKG_j = @(x) sin(acos(R0_j./x)) .* abs(R0_j./(x.^2.*sqrt(1-(R0_j./x).^2))) .* (x>=R0_j & x<=z_max);

        vals_V = fKV_j(z);
        vals_G = fKG_j(z);

        vals_V(~isfinite(vals_V)) = 0;
        vals_G(~isfinite(vals_G)) = 0;

        fKV_mix_raw = fKV_mix_raw + vals_V;
        fKG_mix_raw = fKG_mix_raw + vals_G;

    end

    fKV_mix_raw = fKV_mix_raw / N_R0;
    fKG_mix_raw = fKG_mix_raw / N_R0;

    fKT_z = interp1(K, rho, z, 'linear', 0);
    fKV_z = fKV_mix_raw;
    fKG_z = fKG_mix_raw;

    fKT_z(~isfinite(fKT_z)) = 0;
    fKV_z(~isfinite(fKV_z)) = 0;
    fKG_z(~isfinite(fKG_z)) = 0;

    fKT_z = max(fKT_z,0);
    fKV_z = max(fKV_z,0);
    fKG_z = max(fKG_z,0);

    fKT_z = fKT_z ./ trapz(z,fKT_z);

    F_KT = cumtrapz(z,fKT_z);
    F_KV = cumtrapz(z,fKV_z);
    F_KG = cumtrapz(z,fKG_z);

    F_KT = F_KT ./ max(F_KT);
    F_KV = min(max(F_KV,0),1);
    F_KG = min(max(F_KG,0),1);

    S_KT = 1 - F_KT;
    S_KV = 1 - F_KV;
    S_KG = 1 - F_KG;

    pI_from_V = trapz(z, fKV_z .* S_KT);
    pI_from_G = trapz(z, fKG_z .* S_KT);

    pI_from_V = min(max(pI_from_V,0),1);
    pI_from_G = min(max(pI_from_G,0),1);

    pT_from_V = 1 - pI_from_V;
    pT_from_G = 1 - pI_from_G;

    M = [pI_from_V pI_from_G; pT_from_V pT_from_G];
    P0_vec = [1; 0];

    fZ1 = fKV_z .* S_KT + fKT_z .* S_KV;
    fZ2 = fKG_z .* S_KT + fKT_z .* S_KG;

    fZ1 = max(fZ1,0);
    fZ2 = max(fZ2,0);

    fZ1 = fZ1 ./ trapz(z,fZ1);
    fZ2 = fZ2 ./ trapz(z,fZ2);

    F_Z1 = cumtrapz(z,fZ1);
    F_Z2 = cumtrapz(z,fZ2);

    F_Z1 = F_Z1 ./ max(F_Z1);
    F_Z2 = F_Z2 ./ max(F_Z2);

    mean_pdf = zeros(1,n_segments);
    P_T_markov = zeros(1,n_segments);
    KT_frac_markov = zeros(1,n_segments);

    F_Max_prev = [];
    f_Max_prev = [];

    for seg = 1:n_segments

        P_N = M^(seg-1) * P0_vec;

        P_T_markov(seg) = P_N(1)*pT_from_V + P_N(2)*pT_from_G;
        KT_frac_markov(seg) = mean(P_T_markov(1:seg));

        fS = P_N(1)*fZ1 + P_N(2)*fZ2;
        FS = P_N(1)*F_Z1 + P_N(2)*F_Z2;

        if seg == 1
            F_Max_cur = FS;
            f_Max_cur = fS;
        else
            F_Max_cur = F_Max_prev .* FS;
            f_Max_cur = f_Max_prev .* FS + F_Max_prev .* fS;
        end

        norm_factor = trapz(z,f_Max_cur);
        mean_pdf(seg) = trapz(z,z.*f_Max_cur) / norm_factor;

        F_Max_prev = F_Max_cur;
        f_Max_prev = f_Max_cur;

    end

    w_values(idx) = w;
    y_terminal(idx) = P_T_markov(end);
    y_cumulative(idx) = KT_frac_markov(end);
    mean_pdf_all(idx,:) = mean_pdf;

    pI_from_V_all(idx) = pI_from_V;
    pT_from_V_all(idx) = pT_from_V;
    pI_from_G_all(idx) = pI_from_G;
    pT_from_G_all(idx) = pT_from_G;

    fprintf('K_T=%.4f  R=%.4f  terminal=%.4f  cumulative=%.4f\n', K_T/1e6, w, y_terminal(idx), y_cumulative(idx));
    fprintf('pI_from_V=%.6f  pT_from_V=%.6f  pI_from_G=%.6f  pT_from_G=%.6f\n\n', pI_from_V, pT_from_V, pI_from_G, pT_from_G);

end

log_w = log(w_values);
log_y = log(y_terminal);
coeffs = polyfit(log_w, log_y, 1);
b_fit = coeffs(1);
a_fit = exp(coeffs(2));
w_fit = linspace(min(w_values), max(w_values), 200);
y_fit = a_fit * w_fit.^b_fit;


%%
R_sampling = [1 1.125 1.25 1.5 1.75 2 2.25 2.5 2.75 3];
kappa_sampling = [0.928841 0.8185 0.724181 0.5752 0.467194 0.383824 0.321006 0.272336 0.234893 0.204138];

figure('Color','w','Units','inches','Position',[0.5 0.5 5.5 4]); hold on; box on; grid off;
plot(R_sampling, kappa_sampling, 'ko-', 'LineWidth',2.0, 'MarkerFaceColor','w', 'MarkerSize',6, 'DisplayName','Sampling');
plot(w_values, y_terminal, 'rs-', 'LineWidth',2.0, 'MarkerFaceColor','w', 'MarkerSize',6, 'DisplayName','PDF/Markov terminal');
plot(w_values, y_cumulative, 'bd--', 'LineWidth',1.6, 'MarkerFaceColor','w', 'MarkerSize',5, 'DisplayName','PDF/Markov cumulative');
xlabel('$R$','FontSize',16,'Interpreter','latex');
ylabel('$\kappa$','FontSize',16,'Interpreter','latex');
legend('Location','northeast','FontSize',11,'Interpreter','latex','Box','off');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman');
xlim([0.95 3.05]);
ylim([0 1]);

figure('Color','w','Units','inches','Position',[0.5 0.5 5.5 4]); hold on; box on; grid off;
loglog(w_values, y_terminal, 'ks', 'MarkerFaceColor','k', 'MarkerSize',6, 'HandleVisibility','off');
loglog(w_fit, y_fit, 'r-', 'LineWidth',2, 'DisplayName',sprintf('$\\kappa = %.4fR^{%.4f}$', a_fit, b_fit));
xlabel('$R$','FontSize',16,'Interpreter','latex');
ylabel('$\kappa$','FontSize',16,'Interpreter','latex');
legend('Location','best','FontSize',12,'Interpreter','latex','Box','off');
set(gca,'FontSize',16,'LineWidth',1.2,'TickLabelInterpreter','latex','FontName','Times New Roman','XScale','log','YScale','log');

T_out = table(w_values(:), y_terminal(:), y_cumulative(:), pI_from_V_all(:), pT_from_V_all(:), pI_from_G_all(:), pT_from_G_all(:), mean_pdf_all(:,end)./1e6, ...
    'VariableNames', {'R','kappa_terminal','kappa_cumulative','pI_from_V','pT_from_V','pI_from_G','pT_from_G','K_R_final_MPa'});

disp(T_out);

fprintf('Best fit: kappa = %.4f * R^{%.4f}\n', a_fit, b_fit);