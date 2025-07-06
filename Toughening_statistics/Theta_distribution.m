clear all;
format long;

cd('D:\Data_Storage\Toughening_statistics_2D\GB_N1000_1V4');
Theta_mat1 = load('Theta_mat.mat', 'Theta_mat').Theta_mat;
Theta_mat2 = load('Theta_mat2.mat', 'Theta_mat').Theta_mat;
Theta_mat3 = load('Theta_mat3.mat', 'Theta_mat').Theta_mat;
max_rows = max([size(Theta_mat1, 1), size(Theta_mat2, 1), size(Theta_mat3, 1)]);
Theta_mat1 = [Theta_mat1; zeros(max_rows - size(Theta_mat1, 1), size(Theta_mat1, 2))];
Theta_mat2 = [Theta_mat2; zeros(max_rows - size(Theta_mat2, 1), size(Theta_mat2, 2))];
Theta_mat3 = [Theta_mat3; zeros(max_rows - size(Theta_mat3, 1), size(Theta_mat3, 2))];
Theta_mat_combined = [Theta_mat1, Theta_mat2, Theta_mat3];
Theta_mat_subset = Theta_mat_combined(1:86, :);
Theta_mat_subset = Theta_mat_subset(~isnan(Theta_mat_subset) & Theta_mat_subset ~= 0);
Theta_values = deg2rad(Theta_mat_subset(:));

invalid_values = Theta_values <= 0 | Theta_values > pi/2;
if any(invalid_values)
    Theta_values = Theta_values(~invalid_values);
end

if isempty(Theta_values)
    error('Theta_values is empty after cleaning. No valid data to plot or fit.');
end

cd('C:\Users\sajja\Documents\MATLAB\ASUS\Toughening_statistics_2D');

bin_edges = linspace(0, pi/2, 100);
[counts, edges] = histcounts(Theta_values, bin_edges, 'Normalization', 'pdf');
centers = edges(1:end-1) + diff(edges) / 2;

fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 380, 280]);


histogram(Theta_values, bin_edges, 'Normalization', 'pdf', ...
    'FaceColor', [0.6, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 1.0); hold on;

xlabel('$\Theta_{\rm{I_{V}}}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('$f_{\Theta}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'k');

ft1 = fittype('A ./ (1 + exp(M * ((pi/2 - x).^(1/3) - phi)))','independent','x','coefficients',{'A','M','phi'});
opts1 = fitoptions(ft1);
opts1.StartPoint = [max(counts), -10, 0.5];
opts1.Lower = [0, -Inf, -Inf];
opts1.Upper = [Inf, 0, Inf];
[fit1, gof1] = fit(centers', counts', ft1, opts1);

A = fit1.A;
M = fit1.M;
phi = fit1.phi;
fprintf('Fit1 Coefficients:\nA = %.4f\nM = %.4f\nphi = %.4f\n', A, M, phi);

x_fit = linspace(0, pi/2, 1000);
y_fit = A ./ (1 + exp(M * ((pi/2 - x_fit).^(1/3) - phi)));

% Plot the fitted curve
plot(x_fit, y_fit, '-r', 'LineWidth', 3);

% Add a legend for clarity
% legend('Data Histogram', 'Fitted Curve', 'Location', 'best');


 set(gca, 'FontSize', 16, 'LineWidth', 2, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

x_ticks = linspace(0, pi/2, 5);
xticks(x_ticks);
x_tick_labels = {'$0$', '$\frac{\pi}{8}$', '$\frac{\pi}{4}$', '$\frac{3\pi}{4}$', '$\frac{\pi}{2}$'};
xticklabels(x_tick_labels);
ylim([0 1.15]);


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


xlabel('$\Theta_{\rm{I}}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex');
ylabel('$f_{\Theta}$', 'FontSize', 20, 'FontName', 'Times', 'Interpreter', 'latex', 'Color', 'k');
set(gca, 'TickLabelInterpreter', 'latex');
ylim([0, max(pdf_theta) * 1.1]);
xlim([0, theta_max]);
grid off;
box on;
hold off;
 set(gca, 'FontSize', 16, 'LineWidth', 2, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');


%%
clear all;
theta_max = 90;
x = linspace(0, 90, 1000);
orders = -1.5:.1:1.5;

figure;
hold on;
for n = orders
    theta_poly = theta_max * (1 - (x / 90).^n);
    plot(x, theta_poly, 'LineWidth', 1.5, 'DisplayName', sprintf('n = %d', n));
end

title('\theta(x): Polynomial Decay');
xlabel('x (from 0 to 90)');
ylabel('\theta(x)');

grid on;
hold off;
% figure;
hold on;
for k = orders
    if k == 0
        theta_exp = theta_max * ones(size(x));
    else
        theta_exp = theta_max * (exp(-k * (x / 90)) - exp(-k)) / (1 - exp(-k));
    end
    plot(x, theta_exp, 'LineWidth', 1.5, 'DisplayName', sprintf('k = %d', k));
end
title('\theta(x): Exponential Decay');
xlabel('x (from 0 to 90)');
ylabel('\theta(x)');
legend show;
grid on;
hold off;
%%
theta_max = 1.5;
x = linspace(0, 90, 1000);
orders = 0:4:20;
p_values = linspace(0, 1, 9);  

figure;
for i = 1:length(p_values)
    p = p_values(i);
    subplot(3, 3, i);
    hold on;
    for q = orders
        sigmoid_core = 1 ./ (1 + exp(-q * (x / 90 - p)));
        normalized_sigmoid = (sigmoid_core - sigmoid_core(end)) / (sigmoid_core(1) - sigmoid_core(end));
        theta_sigmoid = theta_max * (normalized_sigmoid);
        plot(x, theta_sigmoid, 'LineWidth', 1.5, 'DisplayName', sprintf('k = %d', q));
    end
    title(sprintf('p = %.1f', p));
    xlabel('x (from 0 to 90)');
    ylabel('\theta(x)');
%     legend show;
    grid on;
    hold off;
end

sgtitle('\theta(x): Sigmoidal Decay for Varying p and q');

%%
% clear all;
% theta_max = 90;
% x = linspace(0, 90, 1000);
% orders = 0:2:10;
% 
% % Taylor expansion at x = 90
% figure;
% hold on;
% for n = orders
%     taylor_approx = zeros(size(x));
%     for k = 1:n
%         if k == 1
%             coef = -1;  % First derivative controls slope
%         else
%             coef = (-1)^(k+1) / factorial(k);  % Higher derivatives control curvature
%         end
%         taylor_approx = taylor_approx + coef * ((x - 90) / 90).^k;  % Normalize by (x/90)
%     end
%     % Scale Taylor series to satisfy boundary condition at x = 0
%     taylor_scaled = theta_max * (taylor_approx - taylor_approx(end)) / (taylor_approx(1) - taylor_approx(end));
%     plot(x, taylor_scaled, 'LineWidth', 1.5, 'DisplayName', sprintf('Order %d', n));
% end
% title('\theta(x): Taylor Series Expansion with Scaling');
% xlabel('x (from 0 to 90)');
% ylabel('\theta(x)');
% legend show;
% grid on;
% hold off;


clear all;
theta_max = 90;
x = linspace(0, 90, 1000);
orders = 1:2:10;
sigma = 20;

figure;
hold on;
for n = orders
    theta_power = theta_max * (1 - (x / 90).^n);
    plot(x, theta_power, 'LineWidth', 1.5, 'DisplayName', sprintf('n = %d', n));
end
title('Power-Law Decay');
xlabel('x');
ylabel('\theta(x)');
legend show;
grid on;
hold off;


% figure;
% theta_trig = theta_max * cos(pi / 2 * x / 90);
% plot(x, theta_trig, 'LineWidth', 1.5);
% title('Trigonometric Decay');
% xlabel('x');
% ylabel('\theta(x)');
% grid on;
% 
% % Logistic Decay
% figure;
% hold on;
% for k = orders
%     theta_logistic = theta_max ./ (1 + exp(k * (x / 90 - 0.5)));
%     plot(x, theta_logistic, 'LineWidth', 1.5, 'DisplayName', sprintf('k = %d', k));
% end
% title('Logistic Decay');
% xlabel('x');
% ylabel('\theta(x)');
% legend show;
% grid on;
% hold off;







