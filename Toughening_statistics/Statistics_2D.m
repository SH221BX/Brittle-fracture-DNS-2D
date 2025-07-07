%% Transgranular

clear all;
theta_max = pi / 4; 
K_range = linspace(1, 1 / cos(theta_max), 1000); 
N_values = [1,3,15,100]; 
num_cols = 1;
num_rows = length(N_values);

fig = figure('Color', 'w');
set(fig, 'Position', [100, 100, 380, 280]);
rng(600); Avg = [];
LowerQuartile = [];
UpperQuartile = [];

for idx = 1:length(N_values)
    N = N_values(idx);
    syms k A
    F = (1 / A) * acos(1 / k); 
    Fn = F^N; 
    fn_sym = diff(Fn, k); 
    
    A_value = theta_max;
    fn = matlabFunction(subs(fn_sym, A, A_value));
    Fn_cdf = matlabFunction(subs(Fn, A, A_value));
      
    a = 1; 
    b = 1 / cos(A_value); 
    
    normalization_constant = integral(fn, a, b);
    fn_normalized = @(k) fn(k) / normalization_constant;
    Fn_cdf_normalized = @(k) Fn_cdf(k) / Fn_cdf(b);
    average_k = integral(@(k) k .* fn_normalized(k), a, b); % E[K]
    E_k2 = integral(@(k) k.^2 .* fn_normalized(k), a, b); % E[K^2]
    variance_k = E_k2 - average_k^2;
    std_k = sqrt(variance_k);
    
    Q1 = fzero(@(k) Fn_cdf_normalized(k) - 0.25, [a, b]); % Find K where CDF = 0.25
    Q3 = fzero(@(k) Fn_cdf_normalized(k) - 0.75, [a, b]); % Find K where CDF = 0.75
    Avg = [Avg; average_k];
    LowerQuartile = [LowerQuartile; Q1];
    UpperQuartile = [UpperQuartile; Q3];
    
%   subplot(num_rows, num_cols, idx);
%   hold on;
    col = {'r','b','c','k'};
    pdf_values = fn_normalized(K_range);
    h_pdf = plot(K_range, pdf_values, '-','Color',col{idx}, 'LineWidth', 3, 'DisplayName', 'PDF'); hold on;
    area(K_range, pdf_values, 'EdgeColor', 'none', 'FaceColor', col{idx}, 'FaceAlpha', 0.25, 'HandleVisibility', 'off');    
    
%     h_EK = xline(average_k, 'k-', 'LineWidth', 2, 'DisplayName', '$E[K]$');
%     h_Q1 = xline(Q1, 'm-', 'LineWidth', 1.5, 'DisplayName', '$Q_1$');
%     h_Q3 = xline(Q3, 'c-', 'LineWidth', 1.5, 'DisplayName', '$Q_3$');  
%   title(sprintf('$N = %d$', N), 'FontSize', 12, 'Interpreter', 'latex');

    xlabel('$\mathrm{K}$', 'FontSize', 14, 'Interpreter', 'latex');
    ylabel('$f_{\mathrm{K_M}}$', 'FontSize', 14, 'Interpreter', 'latex');    
    xlim([a-0.01, b+0.01]);
    ylim([0, 11 * 1.1]);
    grid off;box on;

%     if idx == 2
%     legend([h_pdf, h_EK, h_Q1, h_Q3], {'PDF', '$E[K]$', '$Q_1$', '$Q_3$'}, ...
%            'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
%     hold off; end

legend show;
end
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
grid off;


%%

clear all;
filename = 'R_curve_animation.gif'; 
fig = figure('Color', 'w');

set(fig, 'Position', [200, 300, 380, 280]);
xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\rm {K_{\mathrm{R}}}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
grid off; box on

theta_max = pi / 4; 
N_values = 1:100;  
K_min = 1;
K_max = 1 / cos(theta_max); 

yline(K_max,'LineStyle',':','Color','k','LineWidth',2); hold on;
E_exact = zeros(size(N_values));
E_approx = zeros(size(N_values));
integrand = @(k, N) (acos(1 ./ k).^(N - 1)) ./ sqrt(k.^2 - 1);

for i = 1:length(N_values)
    N = N_values(i);
    E_exact(i) = (N / theta_max^N) * integral(@(k) integrand(k, N), K_min, K_max);   
  % E_approx(i) = K_max - (theta_max * K_max * sqrt(K_max^2 - 1)) / N;

plot(N_values(i), E_exact(i), 'o', 'LineWidth',1,'MarkerFaceColor','r','MarkerEdgeColor','none', 'DisplayName', 'Exact','MarkerSize',4);
hold on;
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
grid off;

ylim([1 1.5]);
frame = getframe(gcf);
img = frame2im(frame);
[imind, cm] = rgb2ind(img, 256); 

% if N == 1 
%     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.25);
% else 
%     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.25);
% end
drawnow;
end

xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\rm{E[K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');


% figure;
% plot(N_values, E_exact, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Exact');
% hold on;
% plot(N_values, E_approx, 'r--x', 'LineWidth', 1.5, 'DisplayName', 'Approximation');
% hold off;

% xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
% ylabel('${E_{\mathrm{max}}[K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
% % legend('Location', 'best');
% grid off; box on;


%%  Intergranular

clc; clear all; close all;
num_samples = 100000;
thetaV = linspace(0,pi/2, num_samples);

A = 1.09530; M = -12.2242; phi = 0.8483;
f_V_original = A ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
integral_f_V = trapz(thetaV, f_V_original);
A = A/ integral_f_V;
f_V = f_V_original / integral_f_V;
    

%%
f_K_unnormalized = @(k) (A ./ (1 + exp(M * ((pi/2 - (acos(1 ./ k))).^(1/3) - phi)))) .* (1 ./ (k .* sqrt(k.^2 - 1)));
normalizing_constant = integral(@(x) f_K_unnormalized(x), 1, Inf);
f_K = @(k) f_K_unnormalized(k) / normalizing_constant;
F_K = @(k) arrayfun(@(kk) integral(@(x) f_K(x), 1.001, kk), k);

kmin = 1.0001;
k_values = linspace(kmin, 1000, num_samples);

f_max = @(k, N) N .* (F_K(k).^(N-1)) .* f_K(k);
F_max = @(k, N) F_K(k).^N;

filename = 'Kmax_animation.gif';
fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);

means = zeros(1, num_samples);

N_values = 1:1:100;

for idx = 1:length(N_values)
    N = N_values(idx);

    pdf_values = f_max(k_values, N);
    cdf_values = F_max(k_values, N);
    norm_const_new = trapz(k_values, pdf_values);
    pdf_values = pdf_values / norm_const_new;
    fn_normalized = @(k) interp1(k_values, pdf_values, k, 'linear', 0);
    average_k = integral(@(k) k .* fn_normalized(k), 1, Inf, 'ArrayValued', true);
    means(N) = average_k;
    col = {'r','b','c','k'};
    % 
    % if N == 1
    %     coll = col{1};
    % elseif N == 3
    %     coll = col{2};
    % elseif N == 15
    %     coll = col{3};
    % else
    %     coll = col{4};
    % end
    % 
    % if N == 1 || N == 3 || N == 15 || N == 100
    % 
    %     plot(k_values, pdf_values, '--','Color',coll, 'LineWidth', 2); hold on;
    % 
    %    area(k_values, pdf_values, 'EdgeColor', 'none', 'FaceColor', coll, 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
    %     % xline(average_k, 'k--', 'LineWidth', 2); hold off;
    %     xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    %    ylabel('$f_{\mathrm{max}}$', 'FontSize', 14, 'Interpreter', 'latex');
    %     % title(['N = ', num2str(N)], 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    % end
    %     xlim([0.8, 7]);
    %    ylim([0,2]);
    % 
    % set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % set(gca, 'ycolor', 'k');
    % drawnow;

end

cd('C:\Users\sajja\Documents\MATLAB\ASUS\Fracture_Model3D\Basic102');
save(means);

%%
fig = figure('Color', 'w');

set(fig, 'Position', [200, 300, 380, 280]);
loglog((1:N), means(1:N), 'o', 'LineWidth',1,'MarkerFaceColor','B','MarkerEdgeColor','none');
xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\rm{E[K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');

set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
grid off;
xlim([0, max(N)]);
ylim([1 5]);  hold on;

cd('C:\Users\sajja\Documents\MATLAB\ASUS\Fracture_Model3D\Basic102');
save(means);

%%

clear all;
theta_max = pi / 2;
theta_range = linspace(0, theta_max, 1000);
pdf_theta = ones(size(theta_range)) / (2 * theta_max);

K_range = linspace(1, 1 / cos(theta_max), 100000);
pdf_K = 1 ./ (theta_max * K_range .* sqrt(K_range.^2 - 1));
cdf_K_values = acos(1 ./ K_range) / theta_max;

fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);

% % subplot(1, 3, 1);
 hold on;
area(theta_range, pdf_theta, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.3);
plot(theta_range, pdf_theta, 'b-', 'LineWidth', 3);

% x_ticks = linspace(-theta_max, theta_max, 5);
% xticks(x_ticks);
% x_tick_labels = {'$-\frac{\pi}{4}$', '$-\frac{\pi}{8}$', '$0$', '$\frac{\pi}{8}$', '$\frac{\pi}{4}$'};
% xticklabels(x_tick_labels);

x_ticks = linspace(0, pi/2, 5);
xticks(x_ticks);
x_tick_labels = {'$0$', '$\frac{\pi}{8}$', '$\frac{\pi}{4}$', '$\frac{3\pi}{4}$', '$\frac{\pi}{2}$'};
xticklabels(x_tick_labels);


xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$f_\theta$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex');
ylim([0, max(pdf_theta) * 1.1]);
xlim([0, theta_max]);
grid off;
box on;
hold on;


%%



% clear all;
cd('D:\Data_Storage\Toughening_statistics_2D\GB_N1000_1V4');

Theta_mat1 = load('Theta_mat.mat', 'Theta_mat').Theta_mat;
Theta_mat2 = load('Theta_mat2.mat', 'Theta_mat').Theta_mat;
Theta_mat3 = load('Theta_mat3.mat', 'Theta_mat').Theta_mat;

KR_mat1 = load('KR_mat.mat', 'KR_mat').KR_mat;
KR_mat2 = load('KR_mat2.mat', 'KR_mat').KR_mat;
KR_mat3 = load('KR_mat3.mat', 'KR_mat').KR_mat;

max_rows = max([size(Theta_mat1, 1), size(Theta_mat2, 1), size(Theta_mat3, 1)]);

Theta_mat1 = [Theta_mat1; zeros(max_rows - size(Theta_mat1, 1), size(Theta_mat1, 2))];
Theta_mat2 = [Theta_mat2; zeros(max_rows - size(Theta_mat2, 1), size(Theta_mat2, 2))];
Theta_mat3 = [Theta_mat3; zeros(max_rows - size(Theta_mat3, 1), size(Theta_mat3, 2))];

KR_mat1 = [KR_mat1; zeros(max_rows - size(KR_mat1, 1), size(KR_mat1, 2))];
KR_mat2 = [KR_mat2; zeros(max_rows - size(KR_mat2, 1), size(KR_mat2, 2))];
KR_mat3 = [KR_mat3; zeros(max_rows - size(KR_mat3, 1), size(KR_mat3, 2))];

Theta_mat_combined = [Theta_mat1, Theta_mat2, Theta_mat3];
KR_mat_combined = [KR_mat1, KR_mat2, KR_mat3];

Theta_mat_subset = Theta_mat_combined(1:50, :);

Theta_mat_subset = Theta_mat_subset(~isnan(Theta_mat_subset));
Theta_values = deg2rad(Theta_mat_subset(:));

% bin_edges = linspace(0, pi/2, 10000);
% [counts, edges] = histcounts(Theta_values, bin_edges, 'Normalization', 'pdf');
% centers = edges(1:end-1) + diff(edges) / 2;
% 
% focus_region = centers >= (70.5 * pi / 180) & centers <= (90 * pi / 180);
% theta_focus = centers(focus_region);
% pdf_focus = counts(focus_region);
% fit_func = @(params, theta) params(1) * (1 - theta / (pi/2)).^params(2);
% 
% initial_guess = [2, 2.57];
% options = optimoptions('lsqcurvefit', 'Display', 'none');
% params = lsqcurvefit(fit_func, initial_guess, theta_focus, pdf_focus, [0, 0], [Inf, Inf], options);
% 
% C = params(1);
% n = params(2);

% figure('Color', 'w');
% histogram(Theta_values, 500, 'Normalization', 'pdf', ...
%     'FaceColor', [0.6, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
%  plot(theta_focus, fit_func(params, theta_focus), 'r-', 'LineWidth', 2);
% xlabel('$\Theta$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
% ylabel('$f_{\Theta}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'k');
% 
% 
% disp(C); disp(n);
% xticks([pi/4, pi/3, 5*pi/12, pi/2]);
% xticklabels({'\pi/4', '\pi/3', '5\pi/12', '\pi/2'});
% grid off;
% set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
% xlim([pi/3.5 pi/2]);
% ylim([0 0.8]);



% figure;
% for i = 20:20
%     N_samples = i;
%     KR_sub = KR_mat_combined(1:N_samples, :);
%     KR_max = max(KR_sub, [], 1);
%     KR_max_clean = KR_max(~isnan(KR_max) & KR_max ~= 0);
%     cd('C:\Users\sajja\Documents\MATLAB\ASUS\Toughening_statistics_2D');
% 
%     histogram(KR_max_clean, 100, 'Normalization', 'pdf'); hold on;
%     xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
%     ylabel('PDF', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'r');
%     title(['N = ', num2str(N_samples)], 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
%     grid off;
%     set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
%     set(gca, 'ycolor', 'k');
% 
% end

% 
% fig = figure('Color','w');
% set(fig,'Position',[100,200,600,150]);
% N_samples_values = [1, 10, 35, 50]; % Define the N values for subplots
% for idx = 1:4
%     N_samples = N_samples_values(idx);
%     KR_sub = KR_mat_combined(1:N_samples, :);
%     KR_max = max(KR_sub, [], 1);
%     KR_max_clean = KR_max(~isnan(KR_max) & KR_max ~= 0);
% 
%     % Calculate statistics
%     mean_KR = mean(KR_max_clean);
%     Q1 = quantile(KR_max_clean, 0.25);
%     Q3 = quantile(KR_max_clean, 0.75);
% 
%     % Create subplot
%     subplot(1, 4, idx);
%     histogram(KR_max_clean, 100, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none'); hold on;
%      xline(mean_KR, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Mean = %.2f', mean_KR));
%     % xline(Q1, 'm--', 'LineWidth', 2, 'DisplayName', sprintf('Q1 = %.2f', Q1));
%     % xline(Q3, 'c--', 'LineWidth', 2, 'DisplayName', sprintf('Q3 = %.2f', Q3));
%     xlim([0 15]);
%     xlabel('K', 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'latex');
%     if idx == 1
%     ylabel('$f_{\mathrm{max}}$', 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'latex');
%     end
%     title(sprintf('$N =%d$ ',N_samples), 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'latex');
%     % legend('show', 'Location', 'best');
%     grid off;
%     set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'FontName', 'Times');
% end
% 
% % Adjust overall figure properties
% set(gcf, 'Color', 'w');


fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);
N_samples_values = 1:1:50;
means = zeros(1, length(N_samples_values));

% Calculate mean values
for idx = 1:length(N_samples_values)
    N_samples = N_samples_values(idx);
    KR_sub = KR_mat_combined(1:N_samples, :);
    KR_max = max(KR_sub, [], 1);
    KR_max_clean = KR_max(~isnan(KR_max) & KR_max ~= 0);
    means(idx) = mean(KR_max_clean);
end

% % Fit the data to a power-law function: means = C * N^exponent
% x_data = N_samples_values(:); % Independent variable (N)
% y_data = means(:); % Dependent variable (means)
% 
% fit_func = @(params, x) params(1) * x.^params(2); % Fit function: C * N^exponent
% initial_guess = [1, 1/2]; % Initial guesses for C and exponent
% options = optimoptions('lsqcurvefit', 'Display', 'off'); % Suppress fitting output
% params = lsqcurvefit(fit_func, initial_guess, x_data, y_data, [0, -Inf], [Inf, Inf], options);
% 
% C = params(1); % Fitted coefficient
% exponent = params(2); % Fitted exponent
% 
% N_theory = linspace(min(N_samples_values), max(N_samples_values), 1000);
% fit_line = fit_func(params, N_theory);
x_data = N_samples_values(:);
y_data = means(:);
fit_func = @(C, x) C * x.^(1/3.512);
initial_guess = 1.00;
options = optimoptions('lsqcurvefit', 'Display', 'off');
C = lsqcurvefit(fit_func, initial_guess, x_data, y_data, 0, Inf, options);
exponent = 1/3.51206;
N_theory = linspace(min(N_samples_values), max(N_samples_values), 1000);
fit_line = fit_func(C, N_theory);

% $ \sim N^{\frac{1}{m+1}}$
% Plot data and fit
loglog(N_samples_values, means, 'ko', 'LineWidth', 2, 'DisplayName', 'DNS', 'MarkerSize', 5,'MarkerFaceColor','k'); hold on;
loglog(N_theory, fit_line, 'c--', 'LineWidth', 3, 'DisplayName', '$\sim N^{\frac{1}{m+1}}$');
xlabel('log(N)', 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('$\rm{log{E[K]}}$', 'FontSize', 14, 'FontName', 'Times', 'Interpreter', 'latex');
grid off;
legend('show', 'Location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontName', 'Times');

% Display fitted parameters in the command window
disp(['Fitted C: ', num2str(C)]);
disp(['Fitted Exponent: ', num2str(exponent)]);






%%




% clc; clear all; close all;

num_samples = 1000000;
thetaV = linspace(0, pi/2, num_samples);
A = 1.0530; M = -12.1242; phi = 0.8483;
A = 1.1099; M = -11.9718; phi = 0.8530;

f_V_original = A ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
 integral_f_V = trapz(thetaV, f_V_original);
A = A / integral_f_V;
f_V = f_V_original / integral_f_V;

f_K_unnormalized = @(k) (A ./ (1 + exp(M * ((pi/2 - acos(1 ./ k)).^(1/3) - phi)))) .* (1 ./ (k .* sqrt(k.^2 - 1)));
normalizing_constant = integral(@(x) f_K_unnormalized(x), 1, Inf);
f_K = @(k) f_K_unnormalized(k) / normalizing_constant;
kmin = 1.001;
k_values = linspace(kmin, 10000, num_samples);
f_K_values = f_K(k_values);
F_K_values = cumtrapz(k_values, f_K_values);
f_max = @(N) N .* (F_K_values.^(N-1)) .* f_K_values;
F_max = @(N) F_K_values.^N;

means = zeros(1, 100);
N_values = 1:50;
for idx = 1:length(N_values)
    N = N_values(idx);
    pdf_values = f_max(N);
    norm_const_new = trapz(k_values, pdf_values);
    pdf_values = pdf_values / norm_const_new;
    average_k = trapz(k_values, k_values .* pdf_values);
    means(N) = average_k;
end
% cd('C:\Users\sajja\Documents\MATLAB\ASUS\Fracture_Model3D\Basic102');
% save('means.mat','means');
% fig = figure('Color', 'w');
% 
% set(fig, 'Position', [200, 300, 380, 280]);
loglog((1:N), means(1:N), '-r', 'LineWidth',3,'MarkerFaceColor','B','MarkerEdgeColor','none','DisplayName','MC');
xlabel('$\rm{log(N)}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\rm{log(E[K])}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');

set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
grid off;
xlim([0, max(N)]);
ylim([1 5]);  hold on;










%%



clear all;

cd('D:\Data_Storage\Toughening_statistics_2D\Com_R1_N1000_1V4\');
Theta_mat1 = load('Theta_mat1.mat', 'Theta_mat').Theta_mat;
Theta_mat2 = load('Theta_mat2.mat', 'Theta_mat').Theta_mat;
KR_mat1 = load('KR_mat1.mat', 'KR_mat').KR_mat;
KR_mat2 = load('KR_mat2.mat', 'KR_mat').KR_mat;
N1 = load('N1.mat', 'N').N;
N2 = load('N2.mat', 'N').N;
max_rows = max([size(Theta_mat1, 1), size(Theta_mat2, 1)]);
Theta_mat1 = [Theta_mat1; zeros(max_rows - size(Theta_mat1, 1), size(Theta_mat1, 2))];
Theta_mat2 = [Theta_mat2; zeros(max_rows - size(Theta_mat2, 1), size(Theta_mat2, 2))];
KR_mat1 = [KR_mat1; zeros(max_rows - size(KR_mat1, 1), size(KR_mat1, 2))];
KR_mat2 = [KR_mat2; zeros(max_rows - size(KR_mat2, 1), size(KR_mat2, 2))];
N_mat1 = [N1; zeros(max_rows - size(KR_mat1, 1), size(KR_mat1, 2))];
N_mat2 = [N2; zeros(max_rows - size(KR_mat2, 1), size(KR_mat2, 2))];
Theta_mat_combined = [Theta_mat1, Theta_mat2];
KR_mat_combined = [KR_mat1, KR_mat2];
N_mat_combined = [N_mat1, N_mat2];
Theta_mat_subset = Theta_mat_combined(1:30, :);
Theta_mat_subset = Theta_mat_subset(~isnan(Theta_mat_subset));
Theta_values = deg2rad(Theta_mat_subset(:));
cd('C:\Users\sajja\Documents\MATLAB\ASUS\Toughening_statistics_2D');


K0 = 1; R = 1; R0 = 1;
A = 1.3230; M = -12.2242; phi = 0.8483;
A = 1.095304825563000; M = -11.970504210371489; phi = 0.848294030387932;
A = 1.3095304825563000; M = -11.970504210371489; phi = 0.8583;

thetaV = linspace(0, pi/2, 10000);
f_V_original = A ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
integral_f_V = trapz(thetaV, f_V_original);
A = A / integral_f_V;
f_V = f_V_original / integral_f_V;



theta_T_max = pi/8;
N_points = 1000; % Reduce points for testing speed

syms theta k

% Define K and Jacobian functions
K_IV = (R0 * K0) / cos(theta);
K_I = (R0 * K0) / cos(theta);
K_T = (R * K0) / cos(theta);

J_IV = (R0 * K0) / (K_IV * sqrt(K_IV^2 - (R0 * K0)^2));
J_I = (R0 * K0) / (K_I * sqrt(K_I^2 - (R0 * K0)^2));
J_T = (R * K0) / (K_T * sqrt(K_T^2 - (R * K0)^2));

f_IV_theta = A / (1 + exp(M * (sqrt(pi/2 - theta) - phi)));
f_I_theta = 1 / (pi / 2);
f_T_theta = 1 / theta_T_max;

f_K_IV = f_IV_theta * J_IV;
f_K_I = f_I_theta * J_I;
f_K_T = f_T_theta * J_T;

f_K_IV_fn = matlabFunction(subs(f_K_IV, theta, acos((R0 * K0) ./ k)));
f_K_I_fn = matlabFunction(subs(f_K_I, theta, acos((R0 * K0) ./ k)));
f_K_T_fn = matlabFunction(subs(f_K_T, theta, acos((R * K0) ./ k)));

normalization_constant_IV = 1;  %integral(@(k) real(f_K_IV_fn(k)), 1, Inf, 'ArrayValued', true);
normalization_constant_I =  1;  %integral(@(k) real(f_K_I_fn(k)), 1, Inf, 'ArrayValued', true);
normalization_constant_T = integral(@(k) real(f_K_T_fn(k)), R * K0, R * K0 / cos(theta_T_max), 'ArrayValued', true);

f_K_IV_normalized = @(k) real(f_K_IV_fn(k)) / normalization_constant_IV;
f_K_I_normalized = @(k) real(f_K_I_fn(k)) / normalization_constant_I;
f_K_T_normalized = @(k) real(f_K_T_fn(k)) / normalization_constant_T;

% Define K range and compute cumulative distributions
K_T_min = 1.001;
K_T_max = R * K0 / cos(theta_T_max);
K_range = linspace(K_T_min, K_T_max, N_points);

F_K_IV_vals = cumtrapz(K_range, f_K_IV_normalized(K_range));
F_K_IV = @(k) interp1(K_range, F_K_IV_vals, k, 'linear', 0);

F_K_I_vals = cumtrapz(K_range, f_K_I_normalized(K_range));
F_K_I = @(k) interp1(K_range, F_K_I_vals, k, 'linear', 0);

F_K_T_vals = cumtrapz(K_range, f_K_T_normalized(K_range));
F_K_T = @(k) interp1(K_range, F_K_T_vals, k, 'linear', 0);

f_min_Iv_T = @(k) arrayfun(@(kk) ...
    (kk >= R * K0 && kk <= R * K0 / cos(theta_T_max)) * ...
    (f_K_IV_normalized(kk) .* (1 - F_K_T(kk)) + f_K_T_normalized(kk) .* (1 - F_K_IV(kk))) + ...
    (kk < R * K0 || kk > R * K0 / cos(theta_T_max)) * f_K_IV_normalized(kk), k);

f_min_I_T = @(k) arrayfun(@(kk) ...
    (kk >= R * K0 && kk <= R * K0 / cos(theta_T_max)) * ...
    (f_K_I_normalized(kk) .* (1 - F_K_T(kk)) + f_K_T_normalized(kk) .* (1 - F_K_I(kk))) + ...
    (kk < R * K0 || kk > R * K0 / cos(theta_T_max)) * f_K_I_normalized(kk), k);

p_V = trapz(K_range, f_K_IV_normalized(K_range) .* (1 - F_K_T(K_range)));
p_G = trapz(K_range, f_K_T_normalized(K_range) .* (1 - F_K_I(K_range)));

f_mix = @(k) (p_V*f_min_Iv_T(k) + (p_G) * f_min_I_T(k));
normalization_constant = trapz(K_range, f_mix(K_range));
f_mix = @(k) f_mix(k) / normalization_constant;

FK_vals = cumtrapz(K_range, f_mix(K_range));
bin_edges = linspace(1, R * K0 / cos(theta_T_max), 1000);
N_factor = 1.001;
filename = 'R12_animation.gif'; 
fig = figure('Color', 'w');

% fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);
means = [];
m_DNS = [];

% N_values = [1,3,15,100];
for N = 1:1:60
    disp(N);
    
% subplot(6,5,N);

%  col = {'r','b','c','k'};
% 
%     if N == 1
%         coll = col{1};
%     elseif N == 3
%         coll = col{2};
%     elseif N == 15
%         coll = col{3};
%     else
%         coll = col{4};
%     end
% 
% if N == 1 || N == 3 || N == 15 || N == 100
    pdf_max = @(k) N * (interp1(K_range, FK_vals, k, 'linear', 0).^(N-1)) .* f_mix(k);
    norm_const = 1;  %trapz(K_range, pdf_max(K_range));
    mean_pdf = trapz(K_range, K_range .* pdf_max(K_range)) / norm_const;
    means = [means;mean_pdf];

%    if N == 100
%         plot(K_range, pdf_max(K_range), '-', 'Color', coll, 'LineWidth', 3, 'DisplayName', sprintf('N = %d', N)); hold on;
%         continue;
%     else
        KR_sub = KR_mat_combined(1:N, :);
        KR_max = max(KR_sub, [], 1);
        KR_max_clean = KR_max(~isnan(KR_max) & KR_max ~= 0);
        m_DNS = [m_DNS;mean(KR_max_clean)];

%         histogram(KR_max_clean, bin_edges, 'Normalization', 'pdf', 'FaceColor', coll, 'EdgeColor', 'w', 'FaceAlpha', 0.2, 'HandleVisibility', 'off'); hold on;
%         ylim([0 70]);
%         set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
%         set(gca, 'ycolor', 'k');
%     end
%     plot(K_range, pdf_max(K_range), '-', 'Color', coll, 'LineWidth', 3, 'DisplayName', sprintf('N = %d', N)); hold on;
% end
% legend('show', 'Interpreter', 'latex');


    % legend('show', 'Interpreter', 'latex');
    % xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
    % ylabel('$f_{\mathrm{max}}$', 'FontSize', 14, 'Interpreter', 'latex');
    % drawnow;

% xlabel('$\rm K$', 'Interpreter', 'latex', 'FontSize', 14);
% 
% ylabel('$f_{\rm {mix}}$', 'Interpreter', 'latex', 'FontSize', 14);
% grid on;
% set(gcf, 'Color', 'w');
% % ylim([0 8]);
% drawnow; pause(0.5);

% frame = getframe(gcf);
% img = frame2im(frame);
% [imind, cm] = rgb2ind(img, 256); 
% 
% if N == 1 
%     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
% else 
%     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
% end

end

fig = figure('Color', 'w');

set(fig, 'Position', [200, 300, 380, 280]);
plot((1:N), means(1:N), '-', 'LineWidth',3,'MarkerFaceColor','B','MarkerEdgeColor','none');hold on;
plot((1:N), m_DNS(1:N), 'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);

xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\rm{E[K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');

set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
grid off;
xlim([0, max(N)]);
% ylim([1 5]);  hold on;


%%


clear all; clc;
A = 1.3095304825563000; M = -11.970504210371489; phi = 0.8583;
thetaV = linspace(0, pi/2, 10000);
f_V_original = A ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
integral_f_V = trapz(thetaV, f_V_original);
A = A / integral_f_V;
f_V = f_V_original / integral_f_V;
mval = M; warning off;


fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);

w_values = []; y_values = [];
R0 = 1;
X = 5000;

for CC = 4:1:4
    Ncp = CC;
for i= 1.001:0.1:3.1

    R = i; n_segments = 40;
    z_min = 1.00001;

    theta_max = (pi/2)*(1/Ncp);
    z_max = R/(cos(theta_max));
    z = linspace(z_min,z_max,X);

    if i == (1.001)
        f_KV = @(x) (A./(1+exp(mval*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).* (x >= R & x <= z_max);
        f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).* (x >= R & x <= z_max);
        f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).* (x >= R & x <= z_max);
        upper1 = z_max;
    else
        f_KV = @(x) (A./(1+exp(mval*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).*(x>=R0 & x<=R);
        f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).*(x>=R & x<=z_max);
        f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<=R);
        upper1 = R;
    end

    F_KV = arrayfun(@(xx) integral(f_KV, R0, xx), z);
    F_KT = arrayfun(@(xx) integral(f_KT, R, min(xx, z_max)), z);
    F_KG = arrayfun(@(xx) integral(f_KG, R0, xx), z);

    F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(zz) f_KT(zz),R,kk),k);
    pV_numerical = integral(@(k) f_KV(k).*F_KT_survival(k),R0,z_max);
    lower1 = R0;    lower2 = R; upper2 = z_max;

    pG_part1 = integral(f_KG,lower1,upper1);
    pG_part2 = integral(@(k) f_KG(k).*F_KT_survival(k),lower2,upper2);

    pG_numerical = pG_part1 + pG_part2;
    M = [pV_numerical pG_numerical; 1-pV_numerical 1-pG_numerical];
    P0 = [1;0];

    fZ1 = f_KV(z).*(1-F_KT) + f_KT(z).*(1-F_KV);
    fZ2 = f_KG(z).*(1-F_KT) + f_KT(z).*(1-F_KG);

    F_Z1 = cumtrapz(z,fZ1); F_Z1 = F_Z1/max(F_Z1);
    F_Z2 = cumtrapz(z,fZ2); F_Z2 = F_Z2/max(F_Z2);

    pT_V = trapz(z,f_KT(z).*(1-F_KV))./trapz(z,fZ1);
    pT_G = trapz(z,f_KT(z).*(1-F_KG))./trapz(z,fZ2);

    f_S = cell(1,n_segments); F_S = cell(1,n_segments);
    F_Max = cell(1,n_segments);
    f_Max = cell(1,n_segments);
    T_count = 0;

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
        T_count = (P_N(1)*pT_V + P_N(2)*pT_G);

        if seg == 40
            loglog(i,(P_N(1)*pT_V + P_N(2)*pT_G),'ks'); hold on;
        end
        % plot(seg,T_count,'ro'); hold on;

    end
    w = i;
    y = (P_N(1)*pT_V + P_N(2)*pT_G);
    w_values = [w_values, w];
    y_values = [y_values, y];
    disp(i); disp(y);

end

log_w = log(w_values);
log_y = log(y_values);

coeffs = polyfit(log_w, log_y, 1);
b = coeffs(1);           % Slope
log_a = coeffs(2);       % Intercept
a = exp(log_a);          % Scaling factor

w_fit = linspace(min(w_values), max(w_values), 100);
y_fit = a * w_fit.^b;

loglog(w_fit, y_fit, '-','Color',rand(1,3), 'LineWidth', 2, 'DisplayName', sprintf('Fit: y = %.4f w^{%.4f}', a, b));
% legend('Data', sprintf('Fit: y = %.4f w^{%.4f}', a, b));
xlabel('$\frac{K_T}{K_{GB}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\frac{N_T}{N_{Total}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
grid off;
hold on;
fprintf('Best fit: y = %.4f * w^{%.4f}\n', a, b);
drawnow;
end
%%

log_w = log(w_values);
log_y = log(y_values);

coeffs = polyfit(log_w, log_y, 1);
b = coeffs(1);           % Slope
log_a = coeffs(2);       % Intercept
a = exp(log_a);          % Scaling factor

w_fit = linspace(min(w_values), max(w_values), 100);
y_fit = a * w_fit.^b;

loglog(w_values,y_values,'ks');hold on;
loglog(w_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Fit: y = %.4f w^{%.4f}', a, b));

grid on; hold on;

fprintf('Best fit: y = %.4f * w^{%.4f}\n', a, b);
xlabel('$\rm{log(R)}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
ylabel('$\rm{log(\kappa)}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');

set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');
annotation('textbox', [0.5 0.7 0.3 0.1], 'String', '$\kappa = \frac{1.00}{R^{3.64}}$', 'Interpreter', 'latex', 'FontSize', 12, 'EdgeColor', 'none');

% figure;
% plot(w_values, y_values, 'o', 'MarkerFaceColor', 'c', 'MarkerSize', 8, 'MarkerEdgeColor', 'k'); hold on;
% plot(w_values, a.* (w_values).^(b), '-c', 'LineWidth', 3);
% 
% xlabel('$\frac{K_T}{K_{GB}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
% ylabel('$\frac{N_T}{N_{Total}}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');



%% VORONOI

% R0 = 1 ; R = 1   DONE DONE DONE DONE

clear; clc; format long;

N = 1e6;  R0 = 1;R = 1; 
n_segments = 1000;

fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);

for CC = 2:1:4

    Ncp = CC;
pV = zeros(1, n_segments);
P_T = zeros(1, n_segments);
Segments = cell(1, n_segments);
MaxSegments = cell(1, n_segments);
thetaV = linspace(0,pi/2, N);

A = 1.09530; M = -12.2242; phi = 0.8483;
f_V_original = A ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
integral_f_V = trapz(thetaV, f_V_original);
A = A/ integral_f_V;
f_V = f_V_original / integral_f_V;
    

fVvals = A ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
fVvals = fVvals ./ trapz(thetaV, fVvals);
FVvals = cumtrapz(thetaV, fVvals);
FVvals = FVvals ./ max(FVvals);
invFV = @(u) interp1(FVvals, thetaV, u, 'linear', pi/2);


theta_max = (pi/2)*(1/Ncp);
R_max = R / cos(theta_max);

z_min = 1.00001;
z_max = R_max;
X = 1000;
z = linspace(z_min, z_max, X);

f_KV = @(x) (A ./ (1 + exp(M * ((pi/2 - acos(R0 ./x)).^(1/3) - phi))) ) .* (R0 ./ x.^2) ./ sqrt(1 - (R0 ./x).^2) .* (x >= R & x <= z_max);
f_KT = @(x) (8 * R) ./ (pi * x.^2 .* sqrt(1 - (R ./ x).^2)) .* (x >= R & x <= z_max);
f_KG = @(x) (2 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (x >= R & x <= z_max);

F_KV = arrayfun(@(xx) integral(f_KV, R0, xx), z);
F_KT = arrayfun(@(xx) integral(f_KT, R, min(xx, z_max)), z);
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
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, z_max);

lower1 = R0; upper1 = R_max; lower2 = R; upper2 = z_max; 

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
%         disp(P_N);
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


    norm_factor = trapz(z,f_Max{seg});
    mean_pdf(seg) = trapz(z,z.*f_Max{seg})/norm_factor;

end


plot(1:n_segments,mean_pdf,'-','LineWidth',2); hold on;
xlabel('N','Interpreter','latex','FontSize',14);
ylabel('$\mathrm{E[K]}$','Interpreter','latex','FontSize',14);
grid off;

end

set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');


%%



clear all; clc;
A = 1.3095304825563000; M = -11.970504210371489; phi = 0.8583;
thetaV = linspace(0, pi/2, 10000);
f_V_original = A ./ (1 + exp(M * ((pi/2 - thetaV).^(1/3) - phi)));
integral_f_V = trapz(thetaV, f_V_original);
A = A / integral_f_V;
f_V = f_V_original / integral_f_V;
mval = M; warning off;


fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);

w_values = []; y_values = [];
R0 = 1;
X = 1000;
Ncp = 4;

for i= 1.001:1:3.1

    R = i; n_segments = 1000;
    z_min = 1.00001;

    theta_max = (pi/2)*(1/Ncp);
    z_max = R/(cos(theta_max));
    z = linspace(z_min,z_max,X);

    if i == (1.001)
        f_KV = @(x) (A./(1+exp(mval*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).* (x >= R & x <= z_max);
        f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).* (x >= R & x <= z_max);
        f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).* (x >= R & x <= z_max);
        upper1 = z_max;
    else
        f_KV = @(x) (A./(1+exp(mval*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).*(x>=R0 & x<=R);
        f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).*(x>=R & x<=z_max);
        f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<=R);
        upper1 = R;
    end


    F_KV = arrayfun(@(xx) integral(f_KV, R0, xx), z);
    F_KT = arrayfun(@(xx) integral(f_KT, R, min(xx, z_max)), z);
    F_KG = arrayfun(@(xx) integral(f_KG, R0, xx), z);


    F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(zz) f_KT(zz),R,kk),k);
    pV_numerical = integral(@(k) f_KV(k).*F_KT_survival(k),R0,z_max);

    lower1 = R0;    lower2 = R; upper2 = z_max;


    pG_part1 = integral(f_KG,lower1,upper1);
    pG_part2 = integral(@(k) f_KG(k).*F_KT_survival(k),lower2,upper2);

    pG_numerical = pG_part1 + pG_part2;
    M = [pV_numerical pG_numerical; 1-pV_numerical 1-pG_numerical];
    P0 = [1;0];

    fZ1 = f_KV(z).*(1-F_KT) + f_KT(z).*(1-F_KV);
    fZ2 = f_KG(z).*(1-F_KT) + f_KT(z).*(1-F_KG);

    F_Z1 = cumtrapz(z,fZ1); F_Z1 = F_Z1/max(F_Z1);
    F_Z2 = cumtrapz(z,fZ2); F_Z2 = F_Z2/max(F_Z2);

    pT_V = trapz(z,f_KT(z).*(1-F_KV))./trapz(z,fZ1);
    pT_G = trapz(z,f_KT(z).*(1-F_KG))./trapz(z,fZ2);

    f_S = cell(1,n_segments); F_S = cell(1,n_segments);
    F_Max = cell(1,n_segments);
    f_Max = cell(1,n_segments);
    T_count = 0;

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

plot(1:n_segments,mean_pdf,'-','LineWidth',2); hold on;
xlabel('N','Interpreter','latex','FontSize',14);
ylabel('$\mathrm{E[K]}$','Interpreter','latex','FontSize',14);
grid off;
end

set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times New Roman');















