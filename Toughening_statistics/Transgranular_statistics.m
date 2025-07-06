%% Animation of PDF and CDF

clear all;
filename = 'N_mean_animation.gif'; 
fig = figure('Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off');

% fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 900, 350]);
Avg = [];
M = [];
LowerQuartile = [];
UpperQuartile = [];

for N = 1:1:100
syms k A

F = (1/A) * acos(1/k);
Fn = F^N;

fn_sym = diff(Fn, k);
A_value = deg2rad(180/8);
fn = matlabFunction(subs(fn_sym, A, A_value));

a = 1;                  
b = 1/cos(A_value);     
k = a:0.000001:b;
i = 0; variance_k = [];

for i = 1:length(k)
    cdy(i) = 1 / deg2rad(22.5) * (acos(1 / k(i)));
    pdy(i) = (1 / deg2rad(22.5)) * (1 / k(i)^2) * (1 / sqrt(1 - (1 / (k(i))^2)));
    NPDY(i) = N * ((cdy(i))^(N - 1)) * (pdy(i));
    CPDY(i) = (cdy(i))^N;
end

normalization_constant = integral(fn, a, b);
fn_normalized = @(k) fn(k) / normalization_constant;
average_k = integral(@(k) k .* fn_normalized(k), a, b);  % E(X)
E_k2 = integral(@(k) k.^2 .* fn_normalized(k), a, b);    % E(X^2)

variance_k = E_k2 - average_k^2;
std_k = sqrt(variance_k);
 

m = b - (b*A_value*(sqrt(b^2 - 1))/N);
M = [M; m];
   
pdf_at_avg_k = fn_normalized(average_k);
cdf_at_avg_k = integral(fn_normalized, a, average_k);

lower_index = find(CPDY >= 0.25, 1, 'first'); % Find index for CDF = 0.25
upper_index = find(CPDY >= 0.75, 1, 'first'); % Find index for CDF = 0.75

lower_quartile = k(lower_index);
upper_quartile = k(upper_index);

Avg = [Avg;average_k];
LowerQuartile = [LowerQuartile; lower_quartile];
UpperQuartile = [UpperQuartile; upper_quartile];

% % Plot mean and quartiles
% plot(N, average_k, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k'); hold on;
% plot(N, lower_quartile, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); hold on;
% plot(N, upper_quartile, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); hold on;
% 
% ylabel('K Values');
% xlabel('N');
% legend('Mean', 'Lower Quartile (25%)', 'Upper Quartile (75%)', 'Location', 'best');
% title('Mean and Quartiles vs. N');
% xlim([0, 105]);
% ylim([min(LowerQuartile)-0.01, max(UpperQuartile)+0.01]);
% grid on;
% set(gca, 'FontSize', 14, 'LineWidth', 1.5);

% plot(N,(cdf_at_avg_k),'o','MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;

% plot(N,average_k,'o','MarkerFaceColor','g','MarkerEdgeColor','k'); hold on;
% yline(b,'Color','r','LineWidth',2);
% ylabel('F(K_{max})');
% xlabel('N');
% xlim([-5 N+10]);
% ylim([1.02 1.1]);

% figure;
% yline(b); % yline(.999);

% plot(N,m,'o','MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
% plot(N,average_k,'o','MarkerFaceColor','g','MarkerEdgeColor','k'); hold on;
% yline(b,'Color','r','LineWidth',2);
% ylabel('Mean');
% xlabel('N');
% xlim([-5 N+10]);
% ylim([1.02 1.1]);

% plot(N,abs(m-average_k),'o','MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
% % plot(N,average_k,'o','MarkerFaceColor','g','MarkerEdgeColor','k'); hold on;
% % yline(b,'Color','r','LineWidth',2);
% ylabel('\Delta');
% xlabel('N');
% xlim([-5 N+10]);
% ylim([1.02 1.1]);

% plot(N,E_k2,'5s');hold on;
% plot(N,variance_k,'ro');hold on;
% plot(N,std_k,'bo');hold on;    
% xline(E_k2,'--r','LineWidth',2); hold on;

subplot(1,2,1);
% Plot PDF and key points
plot(average_k, pdf_at_avg_k, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); 
plot(k, NPDY, '-r', 'LineWidth', 2);
xline(average_k, '-k', 'LineWidth', 2);
xline(lower_quartile, 'm-', 'LineWidth', 2);
xline(upper_quartile, 'c-', 'LineWidth', 2);
hold on;

fill([lower_quartile, upper_quartile, upper_quartile, lower_quartile], ...
     [0, 0, max(NPDY), max(NPDY)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;

xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('PDF', 'FontSize', 20, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
title(['N = ', num2str(N)], 'FontSize', 20, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
xlim([0.99 1.09]);
ylim([-5 max(NPDY) + 20]); 
box on;
set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
set(gca, 'Color', 'none');

subplot(1,2,2);
% Plot CDF and key points
plot(average_k, cdf_at_avg_k, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k'); 
plot(k, CPDY, 'b-', 'LineWidth', 2);
xline(average_k, '-k', 'LineWidth', 2); 
xline(lower_quartile, 'm-', 'LineWidth', 2);
xline(upper_quartile, 'c-', 'LineWidth', 2);
hold on;
% Shade the IQR region in the CDF plot
fill([lower_quartile, upper_quartile, upper_quartile, lower_quartile], ...
     [0, 0, 1, 1], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;
% Configure axes and labels
xlabel('K', 'FontSize', 20, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('CDF', 'FontSize', 20, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
title(['N = ', num2str(N)], 'FontSize', 20, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
xlim([0.99 1.09]);
ylim([-0.05 1.02]); % Ensure CDF range is visible
box on;
set(gca, 'FontSize', 18, 'LineWidth', 2, 'FontName', 'Times');
set(gca, 'Color', 'none');

% frame = getframe(gcf);
% img = frame2im(frame);
% [imind, cm] = rgb2ind(img, 256); 
% 
% if N == 1 
%     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1);
% else 
%     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
% end

end

%% R_curve and IQR

N = 1:length(Avg);

figure;
hold on;
fill([N, fliplr(N)], [LowerQuartile; flipud(UpperQuartile)], ...
     'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
plot(N, Avg, 'k-', 'LineWidth',1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k','MarkerSize',3);

plot(N, LowerQuartile, 'm-', 'LineWidth', 1, 'DisplayName', 'Lower Quartile','MarkerSize',3);
plot(N, UpperQuartile, 'c-', 'LineWidth', 1, 'DisplayName', 'Upper Quartile','MarkerSize',3);
yline(1.0824,'Color','b','LineWidth',2); hold on;
xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('${E_{\mathrm{max}}  [K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
% title('Mean and Quartiles vs. N', 'FontSize', 16, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
% legend({'IQR)', 'Mean', 'Lower Quartile', 'Upper Quartile'}, ...
%        'FontSize', 16, 'Location', 'best');
grid off;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'FontName', 'Times');
hold off;
xlim([-5 110]);

%% Static image of PDF and CDF

clear all;
theta_max = pi / 2; 
theta_range = linspace(-theta_max, theta_max, 1000); % Range for theta
pdf_theta = ones(size(theta_range)) / (2 * theta_max); % Uniform PDF for theta

K_range = linspace(1, 1 / cos(theta_max), 1000); % Range for K
pdf_K = 1 ./ (theta_max * K_range .* sqrt(K_range.^2 - 1)); % PDF for K
cdf_K_values = acos(1 ./ K_range) / theta_max; % CDF for K

% Number of subplots for different N values
N_values = [1, 2, 3, 4, 5, 8, 10, 100, 1000];
num_rows = 3;
num_cols = 3;

fig = figure('Color', 'w');
set(fig, 'Position', [100, 0, 1200, 800]);

for idx = 1:length(N_values)
    N = N_values(idx);

    F = acos(1 ./ K_range) / theta_max;
    f_max = N .* (F.^(N - 1)) .* pdf_K; 
    F_max = F.^N; 

    subplot(num_rows, num_cols, idx);
%     yyaxis left; 
    area(K_range, f_max, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.5); hold on;
    plot(K_range, f_max, 'r-', 'LineWidth', 2);
    set(gca, 'ycolor', 'k');
    ylim([0, max(f_max) * 1.1]);
    xlim([0.99, 1.09]);

%     yyaxis right; 
%     area(K_range, F_max, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.3); hold on;
%     plot(K_range, F_max, 'b--', 'LineWidth', 2);
%     ylim([0, 1.1]);
%     set(gca, 'ycolor', 'b');

    % Title for each subplot
    title(sprintf('$N = %d$', N), 'FontSize', 10, 'Interpreter', 'latex');
    grid off;
    box on;
end

set(gcf, 'Color', 'w');



%% Static image sible cleavage plane of PDF and CDF

clear all;

% Parameters
kmin = 1; % Lower bound for k
kmax = 1000; % Upper bound for k (for plotting purposes)
k_values = linspace(kmin, kmax, 100000); % Discrete k values for plotting
N_values = [1, 10, 50, 100]; % N values
num_rows = 1;
num_cols = 4;

% Functions for f_K and F_K
f_K = @(k) (4/pi) .* (1 ./ (k.^2 .* sqrt(1 - (1 ./ k).^2))); 
F_K = @(k) (2/pi) .* acos(1 ./ k);

fig = figure('Color', 'w');
set(fig, 'Position', [100, 200, 600, 150]);
means = zeros(size(N_values));

% Loop over N values
for idx = 1:length(N_values)
    N = N_values(idx);

    % Define f_max for current N
    f_max = @(k) N .* (F_K(k).^(N-1)) .* f_K(k);

    % Compute normalization constant
    normalization_constant = integral(f_max, kmin, Inf, 'ArrayValued', true);

    % Normalize the PDF
    fn_normalized = @(k) f_max(k) / normalization_constant;

    % Compute the expected value (mean) of k
    average_k = integral(@(k) k .* fn_normalized(k), kmin, Inf, 'ArrayValued', true);
    means(idx) = average_k;

    % Evaluate the normalized PDF for plotting
    pdf_values = fn_normalized(k_values);

    % Plotting
%     subplot(num_rows, num_cols, idx);
%     hold on;

    % Plot shaded area under the PDF
    area((k_values), pdf_values, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.3);

    % Plot the PDF curve
    plot((k_values), pdf_values, 'r-', 'LineWidth', 2); hold on;

    % Add vertical line for mean
% %     xline(average_k, 'k--', 'LineWidth', 2, 'DisplayName', '$E[K]$');

    % Formatting
%     title(sprintf('$N = %d$', N), 'FontSize', 12, 'Interpreter', 'latex');

% text(250, 0.95 * max(pdf_values), ...
%          sprintf('$N = %d$', N), ...
%          'FontSize', 14, 'Interpreter', 'latex', ...
%          'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
%          'BackgroundColor', 'white', 'EdgeColor', 'none');

    xlabel('$\mathrm{K}$', 'FontSize', 12, 'Interpreter', 'latex');
    if idx ==1
    ylabel('$f_{\mathrm{max}}$', 'FontSize', 12, 'Interpreter', 'latex');
    end
    xlim([-50, 500]);
    ylim([0, max(pdf_values) * 1.1]);
    grid off;
    box on;
%     legend({'$E[K]$'}, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
%     hold off;
end

set(gcf, 'Color', 'w');









% 
% figure;
% plot(N_values, means, 'o-', 'LineWidth', 1.5);
% xlabel('N', 'FontSize', 14);
% ylabel('Mean K_{max}', 'FontSize', 14);
% % title('Behavior of Mean K_{max} with K_{max} = \infty', 'FontSize', 16);
% grid on;








%%

clear all;

% Parameters
theta_max = pi / 8; % Example value for theta_max
K_range = linspace(1, 1 / cos(theta_max), 1000); % Range for K
N_values = [1, 20, 50]; % Values of N
num_cols = 1;
num_rows = length(N_values);

fig = figure('Color', 'w');
set(fig, 'Position', [100, 100, 1000, 250]);

% Preallocate arrays for results
Avg = [];
LowerQuartile = [];
UpperQuartile = [];

for idx = 1:length(N_values)
    N = N_values(idx);
    
    % Symbolic expressions for F and f_max
    syms k A
    F = (1 / A) * acos(1 / k); % CDF
    Fn = F^N; % CDF for maximum
    fn_sym = diff(Fn, k); % PDF for maximum
    
    % Substitute actual value of A (theta_max in radians)
    A_value = theta_max;
    fn = matlabFunction(subs(fn_sym, A, A_value));
    Fn_cdf = matlabFunction(subs(Fn, A, A_value)); % CDF for maximum
   
    a = 1; % Lower limit of K
    b = 1 / cos(A_value); % Upper limit of K
    
    % Normalize the PDF
    normalization_constant = integral(fn, a, b);
    fn_normalized = @(k) fn(k) / normalization_constant;
    Fn_cdf_normalized = @(k) Fn_cdf(k) / Fn_cdf(b); % Normalize the CDF
    
    % (E[K]) and variance
    average_k = integral(@(k) k .* fn_normalized(k), a, b); % E[K]
    E_k2 = integral(@(k) k.^2 .* fn_normalized(k), a, b); % E[K^2]
    variance_k = E_k2 - average_k^2;
    std_k = sqrt(variance_k);
    
    % Calculate quartiles (Q1 and Q3) using the normalized CDF
    Q1 = fzero(@(k) Fn_cdf_normalized(k) - 0.25, [a, b]); % Find K where CDF = 0.25
    Q3 = fzero(@(k) Fn_cdf_normalized(k) - 0.75, [a, b]); % Find K where CDF = 0.75
    
    Avg = [Avg; average_k];
    LowerQuartile = [LowerQuartile; Q1];
    UpperQuartile = [UpperQuartile; Q3];
    
    subplot(num_rows, num_cols, idx);
    hold on;
    
    pdf_values = fn_normalized(K_range);
    h_pdf = plot(K_range, pdf_values, 'r-', 'LineWidth', 2, 'DisplayName', 'PDF');
    area(K_range, pdf_values, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
   
    h_EK = xline(average_k, 'k-', 'LineWidth', 2, 'DisplayName', '$E[K]$');
    h_Q1 = xline(Q1, 'm-', 'LineWidth', 1.5, 'DisplayName', '$Q_1$');
    h_Q3 = xline(Q3, 'c-', 'LineWidth', 1.5, 'DisplayName', '$Q_3$');
    
    title(sprintf('$N = %d$', N), 'FontSize', 12, 'Interpreter', 'latex');
    if idx == 3
    xlabel('$\mathrm{K}$', 'FontSize', 14, 'Interpreter', 'latex');
end 
    ylabel('$f_{\mathrm{max}}$', 'FontSize', 14, 'Interpreter', 'latex');
    
    xlim([a-0.01, b+0.01]);
    ylim([0, max(pdf_values) * 1.1]);
    grid off;
    box on;
    if idx == 2
    legend([h_pdf, h_EK, h_Q1, h_Q3], {'PDF', '$E[K]$', '$Q_1$', '$Q_3$'}, ...
           'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
    hold off; end
end

%%
clear all;

filename = 'R_curve_animation.gif'; 
fig = figure('Color', 'w');

set(fig, 'Position', [200, 300, 400, 300]);
xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\rm {K_{\mathrm{R}}}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
grid off; box on

theta_max = pi / 4; 
N_values = 1:100;   
K_min = 1;
K_max = 1 / cos(theta_max); 

yline(K_max,'LineStyle','--','Color','k','LineWidth',2); hold on;
E_exact = zeros(size(N_values));
E_approx = zeros(size(N_values));
integrand = @(k, N) (acos(1 ./ k).^(N - 1)) ./ sqrt(k.^2 - 1);

for i = 1:length(N_values)
    N = N_values(i);
    E_exact(i) = (N / theta_max^N) * integral(@(k) integrand(k, N), K_min, K_max);
   
    % E_approx(i) = K_max - (theta_max * K_max * sqrt(K_max^2 - 1)) / N;

plot(N_values(i), E_exact(i), 'o', 'LineWidth',1,'MarkerFaceColor','r','MarkerEdgeColor','k', 'DisplayName', 'Exact');
hold on;
set(gca, 'FontSize', 15, 'LineWidth', 1.5, 'FontName', 'Times');
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

% figure;
% plot(N_values, E_exact, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Exact');
% hold on;
% plot(N_values, E_approx, 'r--x', 'LineWidth', 1.5, 'DisplayName', 'Approximation');
% hold off;

% xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
% ylabel('${E_{\mathrm{max}}[K]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
% % legend('Location', 'best');
% grid off; box on;

%%


C = [2,3,4,5];
theta_max = pi / C; % Maximum angle in radians
N_values = 5:100;  
K_min = 1;
K_max = 1 / cos(theta_max); % Maximum K value based on theta_max

E_exact = zeros(size(N_values));
E_approx = zeros(size(N_values));
integrand = @(k, N) (acos(1 ./ k).^(N - 1)) ./ sqrt(k.^2 - 1);

for i = 1:length(N_values)
    N = N_values(i);
    E_exact(i) = (N / theta_max^N) * integral(@(k) integrand(k, N), K_min, K_max);
   
    E_approx(i) = K_max - (theta_max * K_max * sqrt(K_max^2 - 1)) / N;
end

figure;
plot(N_values, E_exact, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Exact','MarkerSize',4);
hold on;
xlabel('N', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('${E_{\mathrm{max}}  [k]}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
legend('Location', 'best');
grid off; box on;

%%


C = [4,6,8];
N_values = 1:1000;   
K_min = 1;

E_exact = zeros(length(N_values), length(C));
E_approx = zeros(length(N_values), length(C));
integrand = @(k, N) (acos(1 ./ k).^(N - 1)) ./ sqrt(k.^2 - 1);

figure;
colors = lines(length(C));      
for j = 1:length(C)
    theta_max = pi / C(j);      % Maximum angle in radians for each C
    K_max = 1 / cos(theta_max); % Maximum K value based on theta_max for each C
    for i = 1:length(N_values)
        N = N_values(i);
        E_exact(i, j) = (N / theta_max^N) * integral(@(k) integrand(k, N), K_min, K_max);

        E_approx(i, j) = K_max - (theta_max * K_max * sqrt(K_max^2 - 1)) / N;
    end
    
    semilogx(N_values, E_exact(:, j)/K_max, '-', 'LineWidth', 2.5, 'Color', colors(j, :), 'DisplayName', sprintf('Exact, C = %d', C(j)), 'MarkerSize', 4);
    hold on;
%     yline(K_max, '-', 'LineWidth', 1.5, 'Color', colors(j, :), 'HandleVisibility', 'off');
    % plot(N_values, E_approx(:, j), '--', 'LineWidth', 1.5, 'DisplayName', sprintf('Approximation, C = %d', C(j)), 'MarkerSize', 4);
end

xlabel('log(N)', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
ylabel('$\frac{E_{\mathrm{max}}[k]}{K_{\mathrm{max}}}$', 'FontSize', 18, 'FontName', 'Times', 'FontWeight', 'normal', 'Interpreter', 'latex');
legend('Location', 'best');
grid off; box on;
