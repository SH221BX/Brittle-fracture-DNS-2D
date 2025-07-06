clear all; close all; clc;

N = 1e6;
theta1 = -90 + 180 * rand(N, 1); 
theta2 = -90 + 180 * rand(N, 1); 

abs_theta1 = abs(theta1);
abs_theta2 = abs(theta2);
M = min(abs_theta1, abs_theta2);

m_values = linspace(0, 90, 1000); 
f_M = (90 - m_values) / 4050; 
F_M = 1 - (1 - m_values / 90);
plot(m_values, f_M, 'r-', 'LineWidth', 2);
figure;
subplot(1, 2, 1);
% histogram(M, 'Normalization', 'pdf', 'BinWidth', 1, 'FaceColor', [0.8, 0.8, 1], 'EdgeColor', 'none');
% hold on;
plot(m_values, f_M, 'r-', 'LineWidth', 2);
xlabel('M (degrees)');
ylabel('Probability Density');
title('PDF of M = min(|\theta_1|, |\theta_2|)');

xlim([0, 90]);
ylim([0, max(f_M) * 1.1]);

subplot(1, 2, 2);
[cdf_counts, cdf_edges] = histcounts(M, 'Normalization', 'cdf', 'BinWidth', 1);
cdf_centers = cdf_edges(1:end-1) + diff(cdf_edges)/2;
plot(cdf_centers, cdf_counts, 'b-', 'LineWidth', 1.5);
hold on;
% plot(m_values, F_M, 'r--', 'LineWidth', 2);
xlabel('M (degrees)');
ylabel('Cumulative Probability');

figure;
histogram(theta1, 'Normalization', 'pdf', 'BinWidth', 1, 'FaceColor', [0.6, 0.8, 1], 'EdgeColor', 'none');
hold on;
histogram(theta2, 'Normalization', 'pdf', 'BinWidth', 1, 'FaceColor', [1, 0.6, 0.6], 'EdgeColor', 'none');
xlabel('\theta (degrees)');
ylabel('Probability Density');
title('Uniform Distributions of \theta_1 and \theta_2');
legend('\theta_1', '\theta_2');
grid on;
xlim([-90, 90]);


%%   K distribution for a single segment picking the min from two unifrom rvs

clear all; close all; clc;
k_values = linspace(1, 10, 10000); 

f_K = @(k) (4/pi) .* (1 - (2/pi) .* acos(1 ./ k)) .* (1 ./ (k .* sqrt(k.^2 - 1)));
F_K = @(k) 1 - (1 - (2/pi) .* acos(1 ./ k)).^2;

pdf_values = f_K(k_values);
cdf_values = F_K(k_values);

figure;
subplot(1, 2, 1);
plot(k_values, pdf_values, 'r-', 'LineWidth', 2);
xlabel('K');
ylabel('f_K');
grid off;

subplot(1, 2, 2);
plot(k_values, cdf_values, 'b-', 'LineWidth', 2);
xlabel('K');
ylabel('F_K');
grid off;


%%  K distribution for N segment picking the min from two unifrom rvs
% Static frames

clear all;close all;clc;

k_values = linspace(1, 20, 1000);
f_K = @(k) (4/pi) .* (1 - (2/pi) .* acos(1 ./ k)) .* (1 ./ (k .* sqrt(k.^2 - 1)));
F_K = @(k) 1 - (1 - (2/pi) .* acos(1 ./ k)).^2;
f_max = @(k, N) N .* (F_K(k).^(N-1)) .* f_K(k);
F_max = @(k, N) F_K(k).^N;

figure;
subplot(1, 2, 1); hold on;

for N = 20:20:101
    pdf_values = f_max(k_values, N);
    plot(k_values, pdf_values, 'LineWidth', 2, 'DisplayName', sprintf('N = %d', N));
end

xlabel('K');
ylabel('f_{max}(K)');
title('PDF of K_{max} for Different N');
legend;
grid off;


subplot(1, 2, 2);
hold on;
for N = 20:20:101
    cdf_values = F_max(k_values, N);
    plot(k_values, cdf_values, 'LineWidth', 2, 'DisplayName', sprintf('N = %d', N));
end

xlabel('K');
ylabel('F_{max}(K)');
title('CDF of K_{max} for Different N');
legend;
grid off;



%% PDF_CDF (animated) and R_curve in two figues

clear all;
close all;
clc;
cd(pwd);
k = linspace(1, 1000, 10000);
filename = 'Kmax_animation.gif';

f_K = @(k) (4/pi) .* (1 - (2/pi) .* acos(1 ./ k)) .* (1 ./ (k .* sqrt(k.^2 - 1)));
F_K = @(k) 1 - (1 - (2/pi) .* acos(1 ./ k)).^2;

f_max = @(k, N) N .* (F_K(k).^(N-1)) .* f_K(k);
F_max = @(k, N) F_K(k).^N;

means = zeros(1, 100);

for N = 1:100
    fn = @(k) f_max(k, N);
    normalization_constant = integral(fn, 1, Inf);
    fn_normalized = @(k) fn(k) / normalization_constant;
    average_k = integral(@(k) k .* fn_normalized(k), 1, Inf); 
    means(N) = average_k;

    pdf_values = fn_normalized(k);
    cdf_values = F_max(k, N);

    figure(1);
    clf;

    yyaxis left;
    plot(k, pdf_values, 'r-', 'LineWidth', 2);
    xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('PDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'r');
    title(['N = ', num2str(N)], 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    grid on;
    xlim([1, 50]);
    ylim([0, max(pdf_values) * 1.2]);
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    set(gca, 'ycolor', 'r');

    yyaxis right;
    plot(k, cdf_values, 'b-', 'LineWidth', 2);
    ylabel('CDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'b');
    ylim([0, 1.02]);
    set(gca, 'ycolor', 'b');

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if N == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
  
end

figure(2);
plot(1:100, means, 'k-', 'LineWidth', 2);
xlabel('N', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('Mean of K_{max}', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'tex');
grid on;
set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');

%% %% PDF_CDF and R_curve animation in a single figure

clear all;
close all;
clc;
cd('C:\Users\SH\Documents\Docs\MATLAB\New folder\Fracture_Model2D\Final_2D_Github\Mixed_mode');

k = linspace(1, 10000, 100000);
filename = 'GB_Kmax_animation.gif';
fig = figure('Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');
set(fig, 'Position', [200, 300, 1200, 300]);

f_K = @(k) (4/pi) .* (1 - (2/pi) .* acos(1 ./ k)) .* (1 ./ (k .* sqrt(k.^2 - 1)));
F_K = @(k) 1 - (1 - (2/pi) .* acos(1 ./ k)).^2;

f_max = @(k, N) N .* (F_K(k).^(N-1)) .* f_K(k);
F_max = @(k, N) F_K(k).^N;

means = zeros(1, 100);

for N = 1:1:350
    fn = @(k) f_max(k, N);
    normalization_constant = integral(fn, 1, Inf, 'ArrayValued', true); % Normalize over [1, 100]
    fn_normalized = @(k) fn(k) / normalization_constant;
    average_k = integral(@(k) k .* fn_normalized(k), 1, Inf, 'ArrayValued', true);
    means(N) = average_k;

    pdf_values = fn_normalized(k);
    cdf_values = F_max(k, N);

    clf;

    subplot(1, 3, 1); 
%     yyaxis left;
    plot(k, pdf_values, 'r-', 'LineWidth', 2);
    xline(average_k, 'k--', 'LineWidth', 2); % Mean line
    xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('PDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'r');
   
    grid off;
    xlim([1, 50]);
    ylim([0, max(pdf_values) * 1.2]);
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    set(gca, 'ycolor', 'k');

    subplot(1, 3, 2);
    plot(k, cdf_values, 'b-', 'LineWidth', 2);
     xline(average_k, 'k--', 'LineWidth', 2); % Mean 
      title(['N = ', num2str(N)], 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('CDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'b');
    ylim([0, 1.02]); xlim([1, 50]);
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    set(gca, 'ycolor', 'k');

    subplot(1, 3, 3);
    plot(1:N, means(1:N), 'ko', 'MarkerSize', 3, 'DisplayName', 'Mean'); hold on;
    xlabel('N', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('Mean K_{max}', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'tex');
    grid off;
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    xlim([1, N + 5]);
    ylim([1 20]);
     if N > 1
        N_values = 1:N;
        fit_coeff = polyfit(log(N_values), log(means(1:N)), 1);
        beta = fit_coeff(1);
        N_fit = linspace(1, N, 100);
        mean_fit = exp(polyval(fit_coeff, log(N_fit)));
        plot(N_fit, mean_fit, 'b--', 'LineWidth', 2, 'DisplayName', ['$\sim N^{' num2str(beta, '%.2f') '}$']);
    end
legend('Interpreter', 'latex', 'Location', 'northwest');

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if N == 1 
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0);
    else 
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
    end
  
end

%% K dist of single-pick theta from a sample

clear all;
close all;
clc;
cd('C:\Users\SH\Documents\Docs\MATLAB\New folder\Fracture_Model2D\Final_2D_Github\Mixed_mode');


filename = 'Single_pick_GB_Kmax_animation.gif';
fig = figure('Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');
set(fig, 'Position', [200, 300, 1000, 350]);

kmin = 1; kmax = 1000;
k = linspace(kmin,kmax, 1000);

f_K = @(k) (4/pi) .* (1 ./ (k.^2 .* sqrt(1 - (1 ./ k).^2))); 
F_K = @(k) (2/pi) .* acos(1 ./ k);

f_max = @(k, N) N .* (F_K(k).^(N-1)) .* f_K(k);
F_max = @(k, N) F_K(k).^N;

means = zeros(1, 100);

for N = 1:1:20
    fn = @(k) f_max(k, N);
    normalization_constant = integral(fn,kmin,Inf, 'ArrayValued', true); % Normalize over [1, 100]
    fn_normalized = @(k) fn(k) / normalization_constant;
    average_k = integral(@(k) k .* fn_normalized(k), kmin,Inf, 'ArrayValued', true);
    means(N) = average_k;

    pdf_values = fn_normalized(k);
    cdf_values = F_max(k, N);

    clf;

    subplot(1, 2, 1); 
    yyaxis left;
    plot(k, pdf_values, 'r-', 'LineWidth', 2);
    xline(average_k, 'k--', 'LineWidth', 2); % Mean line
    xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('PDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'r');
    title(['N = ', num2str(N)], 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    grid on;
    xlim([1, 200]);
    ylim([0, max(pdf_values) * 1.2]);
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    set(gca, 'ycolor', 'r');

    yyaxis right;
    plot(k, cdf_values, 'b-', 'LineWidth', 2);
    ylabel('CDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'b');
    ylim([0, 1.02]);
    set(gca, 'ycolor', 'b');

    subplot(1, 2, 2);
    plot(1:N, means(1:N), 'ko', 'MarkerSize', 3, 'DisplayName', 'Mean');
    xlabel('N', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('Mean K_{max}', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'tex');
    grid on;
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    xlim([1, N + 5]);
%     ylim([1 125]);

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if N == 1 
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0);
    else 
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
    end
  
end
%%

clear all;
close all;
clc;

kmin = 1;
kmax_values = [1000, 2000, 5000, 10000];  % Gradually increasing K_max
N = 10;  % Fixed N for demonstration

means = zeros(size(kmax_values));

for idx = 1:length(kmax_values)
    kmax = kmax_values(idx);
    k = linspace(kmin, kmax, 100000);

    f_K = @(k) (4/pi) .* (1 ./ (k.^2 .* sqrt(1 - (1 ./ k).^2))); 
    F_K = @(k) (2/pi) .* acos(1 ./ k);

    f_max = @(k) N .* (F_K(k).^(N-1)) .* f_K(k);

    normalization_constant = integral(f_max, kmin, kmax, 'ArrayValued', true);
    fn_normalized = @(k) f_max(k) / normalization_constant;
    average_k = integral(@(k) k .* fn_normalized(k), kmin, kmax, 'ArrayValued', true);

    means(idx) = average_k;
end

figure;
plot(kmax_values, means, 'o-', 'LineWidth', 1.5);
xlabel('K_{max}', 'FontSize', 14);
ylabel('Mean K_{max}', 'FontSize', 14);
title('Behavior of Mean K_{max} as K_{max} Increases', 'FontSize', 16);
grid on;

%% Single-pick with kmax infinite


clear all;close all;
clc;

kmin = 1; 
N_values = 1:100; 
means = zeros(size(N_values));  

f_K = @(k) (4/pi) .* (1 ./ (k.^2 .* sqrt(1 - (1 ./ k).^2))); 
F_K = @(k) (2/pi) .* acos(1 ./ k);

for idx = 1:length(N_values)
    N = N_values(idx);
    f_max = @(k) N .* (F_K(k).^(N-1)) .* f_K(k);

    % Normalization constant over [1, ∞)
    normalization_constant = integral(f_max, kmin, Inf, 'ArrayValued', true);
    fn_normalized = @(k) f_max(k) / normalization_constant;
    average_k = integral(@(k) k .* fn_normalized(k), kmin, Inf, 'ArrayValued', true);
    means(idx) = average_k;
end

figure;
plot(N_values, means, 'o-', 'LineWidth', 1.5);
xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('${E_{\mathrm{max}}[K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
grid on;


%%  Static frames

clear all;close all;clc;

k_values = linspace(1, 150, 100000);
f_K = @(k) (4/pi) .* (1 ./ (k.^2 .* sqrt(1 - (1 ./ k).^2))); 
F_K = @(k) (2/pi) .* acos(1 ./ k);

% f_K = @(k) 0.0178./(k.^2 .* sqrt(k.^2 - 1));
% F_K = @(k) 0.0178*(180/pi)*sin((pi/180)*acos(1./k));

f_max = @(k, N) N .* (F_K(k).^(N-1)) .* f_K(k);
F_max = @(k, N) F_K(k).^N;

figure;
% subplot(1, 2, 1); 
hold on;
N_f = [1 20 50 100 200];

for i = 1:length(N_f)
    N = N_f(i);
    pdf_values = f_max(k_values, N);
    plot(k_values, pdf_values, 'LineWidth', 2, 'DisplayName', sprintf('N = %d', N));hold on;
end

xlabel('$\mathrm{K}$','Interpreter','latex');
ylabel('$f_{\mathrm{max}}$','Interpreter','latex');
legend;
grid off; box on;
xlim([0 150]);
ylim([0 0.1]);

% subplot(1, 2, 2);
% hold on;
% for N = 20:20:101
%     cdf_values = F_max(k_values, N);
%     plot(k_values, cdf_values, 'LineWidth', 2, 'DisplayName', sprintf('N = %d', N));
% end

% xlabel('K');
% ylabel('F_{max}(K)');
% grid off;


%%

clear all;
close all;
clc;

theta = linspace(-pi/2, pi/2, 1000); 
k = linspace(1, 100, 1000);         

f_theta = @(theta) 1 / pi .* (theta >= -pi/2 & theta <= pi/2); 
f_K = @(k) (4 / pi) .* (1 ./ (k.^2 .* sqrt(1 - (1 ./ k).^2))); 
F_K = @(k) (2 / pi) .* acos(1 ./ k);                   

pdf_theta = f_theta(theta)*ones(size(k,1)); 
pdf_K = f_K(k);            
cdf_K = F_K(k);             

figure;
plot(theta, pdf_theta, 'g-', 'LineWidth', 2);
xlabel('\theta (radians)', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('PDF of \theta', 'FontSize', 14, 'Interpreter', 'tex');
grid on;
ylim([0.3182 0.3184]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

figure;
subplot(1, 2, 1);
plot(k, pdf_K, 'r-', 'LineWidth', 2);
xlabel('K', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('PDF of K', 'FontSize', 14, 'Interpreter', 'latex');

xlim([1, 100]);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

subplot(1, 2, 2);
plot(k, cdf_K, 'b-', 'LineWidth', 2);
xlabel('K', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('CDF of K', 'FontSize', 14, 'Interpreter', 'latex');

xlim([1, 100]);
ylim([0, 1.02]);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);


%%  Animation of Acos(theta) function

clc; clear all; close all;
f_K = @(k) 1./(k.^2.*sqrt(k.^2-1));
F_K = @(k) sqrt(k.^2-1)./k;

f_max = @(k, N) N .* (F_K(k).^(N-1)) .* f_K(k);
F_max = @(k, N) F_K(k).^N;

kmin = 1.001;  
kmax = 50;
k = linspace(kmin, kmax, 10000);

filename = 'AcosTheta_animation.gif';
fig = figure('Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');
set(fig, 'Position', [200, 300, 1000, 350]);

means = zeros(1, 100);

for N = 1:20:300
    fn = @(kk) f_max(kk, N);
       
    normalization_constant = integral(fn, kmin, kmax, 'ArrayValued', true);
    fn_normalized = @(kk) fn(kk) / normalization_constant;
    average_k = integral(@(kk) kk .* fn_normalized(kk), kmin, kmax, 'ArrayValued', true);
    means(N) = average_k;

    pdf_values = round(fn_normalized(k),6);
    cdf_values = F_max(k, N);


    clf;

    subplot(1, 2, 1); 
    yyaxis left;
    plot(k, pdf_values, 'r-', 'LineWidth', 2);
    xline(average_k, 'k--', 'LineWidth', 2);
    xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('PDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'r');
    title(['N = ', num2str(N)], 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    grid on;
    xlim([1, 20]);
    ylim([0, max(pdf_values)*1.2]);
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    set(gca, 'ycolor', 'r');

    yyaxis right;
    plot(k, cdf_values, 'b-', 'LineWidth', 2);
    ylabel('CDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'b');
    ylim([0, 1.02]);
    set(gca, 'ycolor', 'b');

    subplot(1, 2, 2);
    plot(1:N, means(1:N), 'ko', 'MarkerSize', 3);
    xlabel('N', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    ylabel('Mean K_{max}', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'tex');
    grid on;
    set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
    xlim([1, N + 5]);
    ylim([1 12]);

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

%     if N == 1 
%         imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1);
%     else 
%         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
%     end
end







