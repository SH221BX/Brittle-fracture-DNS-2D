clear all;
load("A_T.mat");
data_norm = A_T / 90;
[phat, pci] = betafit(data_norm);
x_deg = linspace(0, 90, 100);
x_norm = x_deg / 90;
beta_pdf_norm = betapdf(x_norm, phat(1), phat(2));
beta_pdf_deg = beta_pdf_norm / 90;

figure;
histogram(A_T, 100, 'Normalization', 'pdf','FaceColor','r','EdgeColor','w'); hold on;
plot(x_deg, beta_pdf_deg, 'k-', 'LineWidth', 2);
xlabel('theta (degree)');
ylabel('PDF');
legend('Histogram', 'Beta fit');
% save('A_T.mat', 'A_T');
% save('fitted_pdf.mat', 'x_deg', 'beta_pdf_deg', 'phat');

% clear all;


load("A_T.mat");

data_rad = A_T * pi / 180;  
data_norm = data_rad / (pi/2); 
[phat, pci] = betafit(data_norm);

x_rad = linspace(0, pi/2, 100);  
x_norm = x_rad / (pi/2);  
beta_pdf_norm = betapdf(x_norm, phat(1), phat(2));
beta_pdf_rad = beta_pdf_norm / (pi/2);  

figure;
histogram(data_rad, 100, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'w'); hold on;
plot(x_rad, beta_pdf_rad, 'k-', 'LineWidth', 2);
xlabel('\theta (radians)');
ylabel('PDF');
legend('Histogram', 'Beta fit');
 

%%
R = 1; alpha = phat(1);

beta_param = phat(2);

k = linspace(R, 100, 100000);
theta_deg = (180/pi) * acos(R ./ k);
f_theta = 1/(90 * beta(alpha, beta_param)) * (theta_deg/90).^(alpha-1) .* (1 - theta_deg/90).^(beta_param-1);

dtheta_deg_dk = (180/pi) * (R ./ (k.^2 .* sqrt(1 - (R./k).^2)));
f_K = f_theta .* abs(dtheta_deg_dk);

figure;
plot(k, f_K, 'b-', 'LineWidth', 2);
xlabel('K');
ylabel('f_K');



%%
R = 1;alpha = phat(1);
beta_param = phat(2);

f_K = @(k) 1./(90*beta(alpha,beta_param)) .* (((180/pi)*acos(R./k)/90).^(alpha-1)) .* ...
    (1 - ((180/pi)*acos(R./k)/90)).^(beta_param-1) .* ...
    abs((180/pi)*(R./(k.^2.*sqrt(1 - (R./k).^2))));

k_min = R + 1e-10; 
norm_const = quadgk(f_K, k_min, Inf, 'RelTol',1e-6, 'AbsTol',1e-12);
f_K_normalized = @(k) f_K(k) / norm_const;

k_max = 1000;  
k_values = linspace(k_min, k_max, 100000);
F_K = cumtrapz(k_values, f_K_normalized(k_values));
F_K = F_K / max(F_K);

N_values = 1:1000;
figure;

for i = 1:length(N_values)
    N = N_values(i);
    f_max_vals = N .* (F_K).^(N-1) .* f_K_normalized(k_values);
    % subplot(ceil(length(N_values)/5), 5, i);
    plot(k_values, f_max_vals,'r-' ,'LineWidth',2 ); hold on;
    area(k_values, f_max_vals, 'EdgeColor', 'none', 'FaceColor', 'c', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
    xlim([1 20]);
  
    title(['$\rm N = \, ' num2str(N) '$'],'Interpreter','latex', 'FontSize',8);
    set(gca, 'FontSize', 8, 'LineWidth', 1.25, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    grid off;
end
% xlabel('$\mathrm{K}$', 'FontSize', 14, 'Interpreter', 'latex');
% ylabel('$f_{\mathrm{K}}$', 'FontSize', 14, 'Interpreter', 'latex');

E_K = zeros(size(N_values));
for i = 1:length(N_values)
    N = N_values(i);
    f_max_vals = N .* (F_K).^(N-1) .* f_K_normalized(k_values);
    E_K(i) = trapz(k_values, k_values .* f_max_vals);
end


figure;
plot(N_values, E_K, 'ro-', 'LineWidth', 2, 'MarkerSize', 3); hold on;
xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\rm{E[K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
grid off;



%%  Transgranular PDF single crystal

clear all; clc; close all;
N = 1e7;
K0 = 1;

n_xy0 = [0; 0; 1];
n_yz0 = [1; 0; 0];
n_zx0 = [0; 1; 0];

minKeq_vals = zeros(N,1);
z_hat = [0; 0; 1];

for i = 1:N
    % alpha = 2*pi*rand;         % rotation about z, uniform in [0, 2pi]
    % beta  = acos(2*rand - 1);  % rotation about y; ensures cos(beta) is uniform in [-1,1]
    % gamma = 2*pi*rand;         % second rotation about z, uniform in [0, 2pi]

    alpha = 2 * pi * rand;  % Rotation about z, uniform in [0, 2pi]
    beta_min = 0;  beta_max = pi/2 ;

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


figure;
% histogram(minKeq_vals, 'Normalization','pdf', 'BinEdges', linspace(1,max(minKeq_vals),200));
% hold on;

R = 1;

theta_max1 = acos(1/sqrt(2));
theta_max2 = acos(1/sqrt(3));

k_min = R ;
k_max = R / cos(theta_max2);

k = linspace(k_min, k_max, 2000);

theta = acos(R./k);

f_theta = zeros(size(theta));
for i = 1:length(k)
    if theta(i) <= theta_max1
        f_theta(i) = 3*sin(theta(i));
    elseif theta(i) <= theta_max2
        y_start = 3*sin(theta_max1);
        y_end   = 0;
        ratio = (theta(i) - theta_max1) / (theta_max2 - theta_max1);
        f_theta(i) = y_start + (y_end - y_start)*sqrt(ratio);
    else
        f_theta(i) = 0;
    end
end

dtheta_dk = 1./(k.*sqrt(k.^2 - R^2));
f_K = f_theta .* abs(dtheta_dk);

plot(k, f_K, 'r-', 'LineWidth',2); hold on;
xlabel('k'); ylabel('f_K(k)');
grid on;







%%  Transgranular PDF Evolution

clear; clc; 
% close all;
R = 1;

theta_max1 = acos(1/sqrt(2));
theta_max2 = acos(1/sqrt(3));
k_min = R +.0000001;
k_max = R/cos(theta_max2);
k = linspace(k_min, k_max, 10000);
theta = acos(R./k);
f_theta = zeros(size(theta));

for i = 1:length(k)
    if theta(i) <= theta_max1
        f_theta(i) = 3*sin(theta(i));
    elseif theta(i) <= theta_max2
        y_start = 3*sin(theta_max1);
        y_end = 0;
        ratio = (theta(i) - theta_max1)/(theta_max2 - theta_max1);
        f_theta(i) = y_start + (y_end - y_start)*sqrt(ratio);
    else
        f_theta(i) = 0;
    end
end

dtheta_dk = 1./(k.*sqrt(k.^2 - R^2));
f_K = f_theta .* abs(dtheta_dk);
F_K = cumtrapz(k, f_K);
F_K = F_K/F_K(end);
N_values = 1:1:10000;


for j = 1:length(N_values)
    N = N_values(j);
    f_max = N*(F_K.^(N-1)).*f_K;
    % subplot(ceil(length(N_values)/5), 5, j);
    % 
    % 
%      plot(k, f_max,'b-' ,'LineWidth',2.5 ); 
    % area(k, f_max, 'EdgeColor', 'none', 'FaceColor', '[0.8 0.8 .8]', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    % 
    % title(['$\rm N = \, ' num2str(N) '$'],'Interpreter','latex', 'FontSize',8);
    % set(gca, 'FontSize', 8, 'LineWidth', 1.25, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % grid off;
    drawnow;
    expVals(j) = trapz(k, k .* f_max);
end

% figure;
plot(N_values, expVals, 'r-', 'LineWidth', 2);
xlabel('N'); ylabel('E[K]');
grid off;


%%

clear all; clc;
R0 = 1; R = 1; a = 1.9896; b = 1.5718;
X = 10000;
theta_max1 = acos(1/sqrt(2));
theta_max2 = acos(1/sqrt(3));
z_min = 1 + 1e-10;
z_maxT = R/cos(theta_max2);
zT = linspace(z_min, z_maxT, X);
theta = acos(R./zT);
f_theta = zeros(size(theta));
for i = 1:length(zT)
    if theta(i) <= theta_max1
        f_theta(i) = 3*sin(theta(i));
    elseif theta(i) <= theta_max2
        y_start = 3*sin(theta_max1);
        y_end = 0;
        ratio = (theta(i) - theta_max1)/(theta_max2 - theta_max1);
        f_theta(i) = y_start + (y_end - y_start)*sqrt(ratio);
    else
        f_theta(i) = 0;
    end
end
dtheta_dk = R ./ (zT.^2 .* sqrt(1 - (R./zT).^2));
f_KT = f_theta .* abs(dtheta_dk);

f_KV = @(x) 1./((pi/2)*beta(a,b)) .* ((acos(R0./x)/(pi/2)).^(a-1)) .* (1 - acos(R0./x)/(pi/2)).^(b-1) .* abs(R0./(x.^2.*sqrt(1 - (R0./x).^2)));
f_KG = @(x) sin(acos(R0./x)) .* abs(R0./(x.^2.*sqrt(1 - (R0./x).^2)));

norm_const = quadgk(f_KV, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
f_K_normalized = @(x) f_KV(x) / norm_const;

figure
subplot(1,3,1)
plot(zT, f_KT, '-k','LineWidth',3)
xlim([z_min z_maxT])
ylim([0 3])
ylabel('$f_{\mathrm{K_T}}$','Interpreter','latex')
xlabel('$\rm{K}$','Interpreter','latex')

zI = linspace(z_min, 10, X);
subplot(1,3,2)
plot(zI, f_KG(zI), '-b','LineWidth',3)
xlim([z_min zI(end)])
% ylim([0 10])
ylabel('$f_{\mathrm{K_I}}$','Interpreter','latex')
xlabel('$\rm{K}$','Interpreter','latex')

subplot(1,3,3)
plot(zI, f_K_normalized(zI), '-r','LineWidth',3)
xlim([z_min zI(end)])
% ylim([0 10])
ylabel('$f_{\mathrm{K_{Iv}}}$','Interpreter','latex')
xlabel('$\rm{K}$','Interpreter','latex')


%%   Markov Chain 

    
clear all; clc;
R0 = 1; R = 1; n_segments = 20;
a = 1.9896; b = 1.5718;


X = 10000;
theta_max1 = acos(1/sqrt(2));   
theta_max2 = acos(1/sqrt(3)); 
z_min = 1 + (1e-10);
z_max = R/cos(theta_max2);
z = linspace(z_min, z_max, X);

theta = acos(R./z);
f_theta = zeros(size(theta));

for i = 1:length(z)
    if theta(i) <= theta_max1
        f_theta(i) = 3*sin(theta(i));
    elseif theta(i) <= theta_max2
        y_start = 3*sin(theta_max1);
        y_end = 0;
        ratio = (theta(i) - theta_max1)/(theta_max2 - theta_max1);
        f_theta(i) = y_start + (y_end - y_start)*sqrt(ratio);
    else
        f_theta(i) = 0;
    end
end

dtheta_dk = (R ./ (z.^2 .* sqrt(1 - (R./z).^2)));
f_KT = f_theta .* abs(dtheta_dk);
F_KT = cumtrapz(z, f_KT);
F_KT = F_KT / F_KT(end);

f_KV = @(x) 1./((pi/2)*beta(a, b)) .*((acos(R0./x)/(pi/2)).^(a-1)) .*(1 - acos(R0./x)/(pi/2)).^(b-1) .*abs(R0./(x.^2.*sqrt(1 - (R0./x).^2)));
f_KG = @(x) sin(acos(R0./x)).*abs(R0./(x.^2.*sqrt(1 - (R0./x).^2)));

norm_const = quadgk(f_KV, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
f_K_normalized = @(x) f_KV(x) / norm_const;
F_KV = cumtrapz(z, f_K_normalized(z));
F_KV = F_KV / max(F_KV);

norm_G = quadgk(f_KG, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
f_KG_normalized = @(x) f_KG(x) / norm_G;
F_KG = cumtrapz(z, f_KG_normalized(z));
F_KG = F_KG / max(F_KG);

% figure;
%  plot(z, f_KT, '-g');  hold on;
%  plot(z,f_K_normalized(z),'-r'); hold on;
%  plot(z,f_KG(z),'-b');  hold on;
%  ylim([0 3]); xlim([0.99 1.74]);



% figure
% subplot(1,3,1)
% plot(z, f_KT, '-g')
% xlim([0.99 1.74])
% ylim([0 3])
% ylabel('$f_{\mathrm{K_T}}$','Interpreter','latex')
% xlabel('$\rm{K}$','Interpreter','latex')
% 
% 
% subplot(1,3,2)
% plot(z, f_KG(z), '-b')
% xlim([0.99 1.74])
% % ylim([0 3])
% ylabel('$f_{\mathrm{K_I}}$','Interpreter','latex')
% xlabel('$\rm{K}$','Interpreter','latex')
% 
% subplot(1,3,3)
% plot(z, f_K_normalized(z), '-r')
% xlim([0.99 1.74])
% % ylim([0 3])
% ylabel('$f_{\mathrm{K_{Iv}}}$','Interpreter','latex')
% xlabel('$\rm{K}$','Interpreter','latex')



 %%

f_KT_func = @(zz) interp1(z, f_KT, zz, 'linear', 0);
F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(f_KT_func, R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, z_max);

lower1 = R0;  upper1 = R;
lower2 = R;   upper2 = R/cos(theta_max2);

pG_part1 = integral(f_KG, lower1, upper2);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

M = [pV_numerical,     pG_numerical;
    1 - pV_numerical, 1 - pG_numerical];

P0 = [1; 0];

fZ1 = f_KV(z) .* (1 - F_KT) + f_KT .* (1 - F_KV);
fZ2 = f_KG(z) .* (1 - F_KT) + f_KT .* (1 - F_KG);

% plot(z,fZ2,"LineWidth",2); hold on;

F_Z1 = cumtrapz(z, fZ1);
F_Z1 = F_Z1 / max(F_Z1);
F_Z2 = cumtrapz(z, fZ2);
F_Z2 = F_Z2 / max(F_Z2);


f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
F_Max = cell(1,n_segments);
f_Max = cell(1,n_segments);

pT_V = trapz(z, f_KT_func(z) .* (1 - F_KV)) ./ trapz(z, fZ1);
pT_G = trapz(z, f_KT_func(z) .* (1 - F_KG)) ./ trapz(z, fZ2);

disp(pT_V); disp(pT_G);


for seg = 1:n_segments

    if seg == 1
        P_N = P0;
    else
        [V0,D0] = eig(M);
        D0(1,1) = D0(1,1)^(seg-1);
        D0(2,2) = D0(2,2)^(seg-1);
        M_pow = V0*D0/V0;
        P_N = M_pow*P0;
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

    if seg == 1
        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end

    % T_count = seg*(P_N(1)*pT_V + P_N(2)*pT_G);
    % if seg == 40
    %     plot(i,(P_N(1)*pT_V + P_N(2)*pT_G),'ro'); hold on;
    % end
    % plot(seg,T_count,'ro'); hold on;

    norm_factor = trapz(z,f_Max{seg});
    mean_pdf(seg) = trapz(z,z.*f_Max{seg})/norm_factor;
    
    % subplot(ceil(n_segments/5),5,seg);
    plot(z,f_Max{seg},'r-','LineWidth',1.5);
    title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    xlim([0.99 1.74]);
    % ylim([0 6]);
    drawnow;
end

% figure;
% plot(1:n_segments,mean_pdf,'k-','LineWidth',2); hold on;
% xlabel('N','Interpreter','latex','FontSize',14);
% ylabel('$\mathrm{E[K_{M}]}$','Interpreter','latex','FontSize',14);
% grid on;





%%  ratio more than 1


clear all;clc;
R0 = 1; R = 2; n_segments = 30;
a = 1.9896; b = 1.5718;

X = 10000;
theta_max1 = acos(1/sqrt(2));   
theta_max2 = acos(1/sqrt(3)); 
z_min = 1 + (1e-10);
z_max = R/cos(theta_max2);
z = linspace(z_min, z_max, X);

theta = acos(R./z);
f_theta = zeros(size(theta));

for i = 1:length(z)
    if theta(i) <= theta_max1
        f_theta(i) = 3*sin(theta(i));
    elseif theta(i) <= theta_max2
        y_start = 3*sin(theta_max1);
        y_end = 0;
        ratio = (theta(i) - theta_max1)/(theta_max2 - theta_max1);
        f_theta(i) = y_start + (y_end - y_start)*sqrt(ratio);
    else
        f_theta(i) = 0;
    end
end

dtheta_dk = real(R ./ (z.^2 .* sqrt(1 - (R./z).^2)));


f_KT = f_theta .* abs(dtheta_dk);
f_KV = @(x) 1./((pi/2)*beta(a, b)) .*((acos(R0./x)/(pi/2)).^(a-1)) .*(1 - acos(R0./x)/(pi/2)).^(b-1) .*abs(R0./(x.^2.*sqrt(1 - (R0./x).^2))) .*(x>=R0 & x<=z_max);
f_KG = @(x) sin(acos(R0./x)).*abs(R0./(x.^2.*sqrt(1 - (R0./x).^2))).*(x>=R0 & x<=z_max);


F_KT = cumtrapz(z, f_KT);
F_KT = F_KT / F_KT(end);

norm_const = quadgk(f_KV, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
f_K_normalized = @(x) f_KV(x) / norm_const;
F_KV = cumtrapz(z, f_K_normalized(z));
F_KV = F_KV / max(F_KV);

norm_G = quadgk(f_KG, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
f_KG_normalized = @(x) f_KG(x) / norm_G;
F_KG = cumtrapz(z, f_KG_normalized(z));
F_KG = F_KG / max(F_KG);

% plot(z, f_KT, '-g'); 
% hold on;
% plot(z,f_K_normalized(z),'-r'); hold on;
% plot(z,f_KG(z),'-b');  hold on;
% ylim([0 3]); xlim([0.99 1.74]);

f_KT_func = @(zz) interp1(z, f_KT, zz, 'linear', 0);
F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(f_KT_func, R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), z_min, z_max);

lower1 = R0;  upper1 = R;
lower2 = R;   upper2 = R/cos(theta_max2);

pG_part1 = integral(f_KG, lower1, upper2);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

M = [pV_numerical,     pG_numerical;
    1 - pV_numerical, 1 - pG_numerical];

P0 = [1; 0];

fZ1 = f_KV(z) .* (1 - F_KT) + f_KT .* (1 - F_KV);
fZ2 = f_KG(z) .* (1 - F_KT) + f_KT .* (1 - F_KG);

F_Z1 = cumtrapz(z, fZ1);
F_Z1 = F_Z1 / max(F_Z1);
F_Z2 = cumtrapz(z, fZ2);
F_Z2 = F_Z2 / max(F_Z2);


f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
F_Max = cell(1,n_segments);
f_Max = cell(1,n_segments);


pT_V = trapz(z, f_KT_func(z) .* (1 - F_KV)) ./ trapz(z, fZ1);
pT_G = trapz(z, f_KT_func(z) .* (1 - F_KG)) ./ trapz(z, fZ2);

disp(pT_V); disp(pT_G);

for seg = 1:n_segments

    if seg == 1
        P_N = P0;
    else
        [V0,D0] = eig(M);
        D0(1,1) = D0(1,1)^(seg-1);
        D0(2,2) = D0(2,2)^(seg-1);
        M_pow = V0*D0/V0;
        P_N = M_pow*P0;
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

    if seg == 1
        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end

    T_count = seg*(P_N(1)*pT_V + P_N(2)*pT_G);
    % if seg == length(n_segments)
    %     plot(seg,(P_N(1)*pT_V + P_N(2)*pT_G),'ro'); hold on;
    % end
    plot(seg,T_count,'ro'); hold on;

    norm_factor = trapz(z,f_Max{seg});
    mean_pdf(seg) = trapz(z,z.*f_Max{seg})/norm_factor;
    
    % subplot(ceil(n_segments/5),5,seg);
    % plot(z,f_Max{seg},'r-','LineWidth',1.5);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % % xlim([0.99 1.74]);
    % % ylim([0 6]);
    % drawnow;
end

figure;
plot(1:n_segments,mean_pdf,'c-','LineWidth',2);
xlabel('N','Interpreter','latex','FontSize',14);
ylabel('$\mathrm{E[K_{M}]}$','Interpreter','latex','FontSize',14);
grid on;

%%



clear all;
clc;


for w = 1:0.1:3

R0 = 1; R = w; n_segments = 100;
a = 1.9896; b = 1.5718;


X = 10000;
theta_max1 = acos(1/sqrt(2));   
theta_max2 = acos(1/sqrt(3)); 
z_min = 1 + (1e-10);
z_max = R/cos(theta_max2);
z = linspace(z_min, z_max, X);

theta = acos(R./z);
f_theta = zeros(size(theta));

for i = 1:length(z)
    if theta(i) <= theta_max1
        f_theta(i) = 3*sin(theta(i));
    elseif theta(i) <= theta_max2
        y_start = 3*sin(theta_max1);
        y_end = 0;
        ratio = (theta(i) - theta_max1)/(theta_max2 - theta_max1);
        f_theta(i) = y_start + (y_end - y_start)*sqrt(ratio);
    else
        f_theta(i) = 0;
    end
end

dtheta_dk = real(R ./ (z.^2 .* sqrt(1 - (R./z).^2)));


f_KT = f_theta .* abs(dtheta_dk);
f_KV = @(x) 1./((pi/2)*beta(a, b)) .*((acos(R0./x)/(pi/2)).^(a-1)) .*(1 - acos(R0./x)/(pi/2)).^(b-1) .*abs(R0./(x.^2.*sqrt(1 - (R0./x).^2))) .*(x>=R0 & x<=z_max);
f_KG = @(x) sin(acos(R0./x)).*abs(R0./(x.^2.*sqrt(1 - (R0./x).^2))).*(x>=R0 & x<=z_max);


F_KT = cumtrapz(z, f_KT);
F_KT = F_KT / F_KT(end);

norm_const = quadgk(f_KV, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
f_K_normalized = @(x) f_KV(x) / norm_const;
F_KV = cumtrapz(z, f_K_normalized(z));
F_KV = F_KV / max(F_KV);

norm_G = quadgk(f_KG, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
f_KG_normalized = @(x) f_KG(x) / norm_G;
F_KG = cumtrapz(z, f_KG_normalized(z));
F_KG = F_KG / max(F_KG);

% plot(z, f_KT, '-g'); 
% hold on;
% plot(z,f_K_normalized(z),'-r'); hold on;
% plot(z,f_KG(z),'-b');  hold on;
% ylim([0 3]); xlim([0.99 1.74]);

f_KT_func = @(zz) interp1(z, f_KT, zz, 'linear', 0);
F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(f_KT_func, R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), z_min, z_max);

lower1 = R0;  upper1 = R;
lower2 = R;   upper2 = R/cos(theta_max2);

pG_part1 = integral(f_KG, lower1, upper2);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

M = [pV_numerical,     pG_numerical;
    1 - pV_numerical, 1 - pG_numerical];

P0 = [1; 0];

fZ1 = f_KV(z) .* (1 - F_KT) + f_KT .* (1 - F_KV);
fZ2 = f_KG(z) .* (1 - F_KT) + f_KT .* (1 - F_KG);

F_Z1 = cumtrapz(z, fZ1);
F_Z1 = F_Z1 / max(F_Z1);
F_Z2 = cumtrapz(z, fZ2);
F_Z2 = F_Z2 / max(F_Z2);


f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
F_Max = cell(1,n_segments);
f_Max = cell(1,n_segments);


pT_V = trapz(z, f_KT_func(z) .* (1 - F_KV)) ./ trapz(z, fZ1);
pT_G = trapz(z, f_KT_func(z) .* (1 - F_KG)) ./ trapz(z, fZ2);

% disp(pT_V); disp(pT_G);

colo = 1-rand(1,3);
for seg = 1:n_segments

    if seg == 1
        P_N = P0;
    else
        [V0,D0] = eig(M);
        D0(1,1) = D0(1,1)^(seg-1);
        D0(2,2) = D0(2,2)^(seg-1);
        M_pow = V0*D0/V0;
        P_N = M_pow*P0;
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

    if seg == 1
        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end

    T_count = seg*(P_N(1)*pT_V + P_N(2)*pT_G);
    % if seg == length(n_segments)
    %     plot(seg,(P_N(1)*pT_V + P_N(2)*pT_G),'ro'); hold on;
    % end
    % plot(seg,T_count,'s','Color',colo); hold on;

    norm_factor = trapz(z,f_Max{seg});
    mean_pdf(seg) = trapz(z,z.*f_Max{seg})/norm_factor;
    
    % subplot(ceil(n_segments/5),5,seg);
    % plot(z,f_Max{seg},'r-','LineWidth',1.5);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % % xlim([0.99 1.74]);
    % % ylim([0 6]);
    % drawnow;
end

 loglog(w,(P_N(1)*pT_V + P_N(2)*pT_G),'-ks'); hold on;
end


%%


clear all; clc;

fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);

w_values = []; y_values = [];

for w = 1:0.1:3
    R0 = 1; R = w; n_segments = 100;
    a = 1.9896; b = 1.5718;

    X = 10000;
    theta_max1 = acos(1/sqrt(2));   
    theta_max2 = acos(1/sqrt(3)); 
    z_min = 1 + (1e-10);
    z_max = R/cos(theta_max2);
    z = linspace(z_min, z_max, X);

    theta = acos(R./z);
    f_theta = zeros(size(theta));

    for i = 1:length(z)
        if theta(i) <= theta_max1
            f_theta(i) = 3*sin(theta(i));
        elseif theta(i) <= theta_max2
            y_start = 3*sin(theta_max1);
            y_end = 0;
            ratio = (theta(i) - theta_max1)/(theta_max2 - theta_max1);
            f_theta(i) = y_start + (y_end - y_start)*sqrt(ratio);
        else
            f_theta(i) = 0;
        end
    end

    dtheta_dk = real(R ./ (z.^2 .* sqrt(1 - (R./z).^2)));

    f_KT = f_theta .* abs(dtheta_dk);
    f_KV = @(x) 1./((pi/2)*beta(a, b)) .* ((acos(R0./x)/(pi/2)).^(a-1)) .* (1 - acos(R0./x)/(pi/2)).^(b-1) .* abs(R0./(x.^2.*sqrt(1 - (R0./x).^2))) .* (x>=R0 & x<=z_max);
    f_KG = @(x) sin(acos(R0./x)).*abs(R0./(x.^2.*sqrt(1 - (R0./x).^2))).*(x>=R0 & x<=z_max);

    F_KT = cumtrapz(z, f_KT);
    F_KT = F_KT / F_KT(end);

    norm_const = quadgk(f_KV, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
    f_K_normalized = @(x) f_KV(x) / norm_const;
    F_KV = cumtrapz(z, f_K_normalized(z));
    F_KV = F_KV / max(F_KV);

    norm_G = quadgk(f_KG, z_min, Inf, 'RelTol',1e-8, 'AbsTol',1e-12);
    f_KG_normalized = @(x) f_KG(x) / norm_G;
    F_KG = cumtrapz(z, f_KG_normalized(z));
    F_KG = F_KG / max(F_KG);

    f_KT_func = @(zz) interp1(z, f_KT, zz, 'linear', 0);
    F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(f_KT_func, R, kk), k);
    pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), z_min, z_max);

    lower1 = R0;  upper1 = R;
    lower2 = R;   upper2 = R/cos(theta_max2);

    pG_part1 = integral(f_KG, lower1, upper2);
    pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
    pG_numerical = pG_part1 + pG_part2;

    M = [pV_numerical,     pG_numerical;
         1 - pV_numerical, 1 - pG_numerical];

    P0 = [1; 0];

    fZ1 = f_KV(z) .* (1 - F_KT) + f_KT .* (1 - F_KV);
    fZ2 = f_KG(z) .* (1 - F_KT) + f_KT .* (1 - F_KG);

    F_Z1 = cumtrapz(z, fZ1);
    F_Z1 = F_Z1 / max(F_Z1);
    F_Z2 = cumtrapz(z, fZ2);
    F_Z2 = F_Z2 / max(F_Z2);

    f_S = cell(1,n_segments);
    F_S = cell(1,n_segments);
    F_Max = cell(1,n_segments);
    f_Max = cell(1,n_segments);

    pT_V = trapz(z, f_KT_func(z) .* (1 - F_KV)) ./ trapz(z, fZ1);
    pT_G = trapz(z, f_KT_func(z) .* (1 - F_KG)) ./ trapz(z, fZ2);

    for seg = 1:n_segments
        if seg == 1
            P_N = P0;
        else
            [V0,D0] = eig(M);
            D0(1,1) = D0(1,1)^(seg-1);
            D0(2,2) = D0(2,2)^(seg-1);
            M_pow = V0*D0/V0;
            P_N = M_pow*P0;
        end

        f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
        F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

        if seg == 1
            F_Max{seg} = F_S{seg};
            f_Max{seg} = f_S{seg};
        else
            F_Max{seg} = F_Max{seg-1}.*F_S{seg};
            f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
        end

        norm_factor = trapz(z,f_Max{seg});
        mean_pdf(seg) = trapz(z,z.*f_Max{seg})/norm_factor;
    end

    y = P_N(1)*pT_V + P_N(2)*pT_G;
    w_values = [w_values, w];
    y_values = [y_values, y];
    loglog(w, y, '-ks'); hold on;
end

% xlabel('N','Interpreter','latex','FontSize',14);
% ylabel('$\mathrm{E[K]}$','Interpreter','latex','FontSize',14);
% set(gca, 'FontSize', 16, 'LineWidth', 2, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

log_w = log(w_values);
log_y = log(y_values);

coeffs = polyfit(log_w, log_y, 1);
b = coeffs(1);           % Slope
log_a = coeffs(2);       % Intercept
a = exp(log_a);          % Scaling factor

w_fit = linspace(min(w_values), max(w_values), 100);
y_fit = a * w_fit.^b;

plot(w_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Fit: y = %.4f w^{%.4f}', a, b));
legend('Data', sprintf('Fit: y = %.4f w^{%.4f}', a, b));
xlabel('w');
ylabel('P_N(1)*pT_V + P_N(2)*pT_G');
title('Log-Log Plot with Power-Law Fit');
grid on;

fprintf('Best fit: y = %.4f * w^{%.4f}\n', a, b);
%%
w = 1:0.1:3;
figure;
plot(w_values, y_values, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'MarkerEdgeColor', 'k'); hold on;
plot(w, 0.7295 .* (w).^(-1.9888), '-c', 'LineWidth', 3);
% yline(0.441);hold on;
% yline(0.469); hold on;
xlabel('$\rm{R = \frac{K_T}{K_{GB}}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\rm{\kappa =\frac{N_T}{N_{Total}}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');

%%

w = 1:0.1:3;
fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);
curve = 0.7295 .* (w).^(-1.9887);
loglog(w, curve, '-c', 'LineWidth', 3); hold on;
loglog(w_values, y_values, 'o', 'MarkerFaceColor', 'c', 'MarkerSize', 5, 'MarkerEdgeColor', 'k');
hold on;

y1 = 0.441;
y2 = 0.469;
ym = 0.459;

% yline(y1, 'k--', 'LineWidth',1);
% yline(y2, 'k--', 'LineWidth',1);
% yline(ym, 'k-', 'LineWidth',2);

w_int1 = (0.7295 / y1)^(1/1.9887);
w_int2 = (0.7295 / y2)^(1/1.9887);
w_m = (0.7295 / ym)^(1/1.9887);

% xline(w_int1, 'r--'); xline(w_int2, 'r--'); xline(w_m, 'r-','LineWidth',2);

% xlabel('$\rm{R = \frac{K_T}{K_{GB}}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
% ylabel('$\rm{\kappa = \frac{N_T}{N_{Total}}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');

xlabel('$\rm{log(R)}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\rm{log(\kappa)}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');

set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
annotation('textbox', [0.5 0.7 0.3 0.1], 'String', '$\kappa = \frac{0.7295}{R^2}$', 'Interpreter', 'latex', 'FontSize', 12, 'EdgeColor', 'none');


%%   TaC Experimental Indentation

clear all;  clc;
N_T = [74 98 247 170 146 172];
N_G = [91 115 283 192 185 205];
crack_start = [57 50 99 58 45 59];
Indent = [10 10 20 10 10 10];
Load = [2 4 6 8 10 12];

figure;

R_T = N_T./(N_T+N_G);

Rmin = [0.10 0.235 0.375 0.333 0.291 0.380];
Rmax = [0.61 0.607 0.580 0.585 0.666 0.608];
lower = mean(R_T) - Rmin;
upper = Rmax - mean(R_T);

b = plot(Load, R_T,'ks','MarkerFaceColor','r','MarkerSize',10);hold on;

xtips = b.XData;
errorbar(xtips, mean(R_T), lower, upper, 'k.', ...
    'LineWidth', 1, 'CapSize', 12);

yline(mean(R_T),'Color','r','LineWidth',2,'LineStyle','--'); hold on;

xlabel('$\rm{Load(N)}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\rm{Ratio}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');

ylim([0 0.8]);
xlim([1 13]);


%%
Nx = (N_T + N_G) ./ Indent;
Nx_min = [13 15 16 21 22 27];            
Nx_max = [19 32 35 52 50 52];   

Nt_min = [2 4 6 7 7 12];
Nt_max = [10 17 18 27 22 28];

Ng_min = [5 8 10 14 11 15];
Ng_max =[17 16 20 31 30 32];

lower_err = Nx - Nx_min;
upper_err = Nx_max - Nx;

figure;
b = bar(Load, Nx, 0.6, 'FaceColor', 'flat','EdgeColor', 'k','LineWidth',2,'FaceAlpha',0.6);hold on;
colors = lines(length(Nx));  
b.CData = colors;

xtips = b.XData;
errorbar(xtips, Nx, lower_err, upper_err, 'k.', ...
    'LineWidth', 1.5, 'CapSize', 12);
xlabel('Load (N)', 'FontName','Times New Roman');
ylabel('Average N per sample', 'FontName','Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');

figure;
Ntx = N_T./Indent;
lower_t = Ntx - Nt_min;
upper_t = Nt_max - Ntx;

b = bar(Load, N_T./Indent, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'r','LineWidth',2,'FaceAlpha',0.6);hold on;
colors = lines(length(Nx));  
b.CData = colors;

xtips = b.XData;
errorbar(xtips, Ntx, lower_t, upper_t, 'r.', ...
    'LineWidth', 1.5, 'CapSize', 12);

xlabel('Load (N)', 'FontName','Times New Roman');
ylabel('Average N_T per sample', 'FontName','Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');


figure;
Ngx = N_G./Indent;
lower_g = Ngx - Ng_min;
upper_g = Ng_max - Ngx;

b = bar(Load, N_G./Indent, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'b','LineWidth',2,'FaceAlpha',0.6);hold on;
colors = lines(length(Nx));  
b.CData = colors;
xtips = b.XData;
errorbar(xtips, Ngx, lower_g, upper_g, 'b.', ...
    'LineWidth', 1.5, 'CapSize', 12);
xlabel('Load (N)', 'FontName','Times New Roman');
ylabel('Average N_G per sample', 'FontName','Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');





%%

clear all; clc;

N_T = [74 98 247 170 146 172];
N_G = [91 115 283 192 185 205];
crack_start = [57 50 99 58 45 59];
Indent = [10 10 20 10 10 10];
Load = [2 4 6 8 10 12];

Nx = (N_T + N_G) ./ Indent;
Nx_min = [13 15 16 21 22 27];
Nx_max = [19 32 35 52 50 52];

Nt_min = [2 4 6 7 7 12];
Nt_max = [10 17 18 27 22 28];

Ng_min = [5 8 10 14 11 15];
Ng_max = [17 16 20 31 30 32];

lower_err = Nx - Nx_min;
upper_err = Nx_max - Nx;

figure;

% 1st subplot: Nx
subplot(2,2,1);
b = bar(Load, Nx, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'k','LineWidth',2,'FaceAlpha',0.6); hold on;
b.CData = lines(length(Nx));
errorbar(b.XData, Nx, lower_err, upper_err, 'k.', 'LineWidth', 1.5, 'CapSize', 12);
xlabel('Load (N)', 'FontName','Times New Roman');
ylabel('N per indent', 'FontName','Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
title('a', 'FontName','Times New Roman');

% 2nd subplot: N_T
subplot(2,2,2);
Ntx = N_T ./ Indent;
lower_t = Ntx - Nt_min;
upper_t = Nt_max - Ntx;
b = bar(Load, Ntx, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'r','LineWidth',2,'FaceAlpha',0.6); hold on;
b.CData = lines(length(Nx));
errorbar(b.XData, Ntx, lower_t, upper_t, 'r.', 'LineWidth', 1.5, 'CapSize', 12);
xlabel('Load (N)', 'FontName','Times New Roman');
ylabel('N_T per indent', 'FontName','Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
title('b', 'FontName','Times New Roman');

% 3rd subplot: N_G
subplot(2,2,3);
Ngx = N_G ./ Indent;
lower_g = Ngx - Ng_min;
upper_g = Ng_max - Ngx;
b = bar(Load, Ngx, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'b','LineWidth',2,'FaceAlpha',0.6); hold on;
b.CData = lines(length(Nx));
errorbar(b.XData, Ngx, lower_g, upper_g, 'b.', 'LineWidth', 1.5, 'CapSize', 12);
xlabel('Load (N)', 'FontName','Times New Roman');
ylabel('N_G per indent', 'FontName','Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
title('c', 'FontName','Times New Roman');

% 4th subplot: Ratio plot
subplot(2,2,4);

R_T = N_T./(N_T+N_G);

Rmin = [0.10 0.235 0.375 0.333 0.291 0.380];
Rmax = [0.61 0.607 0.580 0.585 0.666 0.608];
lower = mean(R_T) - Rmin;
upper = Rmax - mean(R_T);

b = plot(Load, R_T,'ks','MarkerFaceColor','r','MarkerSize',8);hold on;

xtips = b.XData;
errorbar(xtips, mean(R_T), lower, upper, 'k.', ...
    'LineWidth', 1, 'CapSize', 12);

yline(mean(R_T),'Color','r','LineWidth',2,'LineStyle','--'); hold on;

xlabel('$\rm{Load(N)}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\rm{\frac{N_{\rm{T}}}{N}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
title('d', 'FontName','Times New Roman');
ylim([0 0.8]);
xlim([1 13]);