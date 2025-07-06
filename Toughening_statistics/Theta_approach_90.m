clear;
close all;
clc;

n_values = [0.5, 1, 2]; 
theta_max = pi/2; 
num_points = 1000;
theta = linspace(0, theta_max, num_points); 

figure;
hold on;
grid on;
xlabel('$\theta$ (radians)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$f_\theta(\theta)$', 'Interpreter', 'latex', 'FontSize', 12);
legend_entries = cell(length(n_values),1);

colors = lines(length(n_values));

for i = 1:length(n_values)
    n = n_values(i);    
    C = ((n + 1) * 2^(n + 1)) / (pi^(n + 1));  
    f_theta = C * ( (pi/2) - theta ).^n;
    f_theta(theta > pi/2) = 0;
    plot(theta, f_theta, 'LineWidth', 2, 'Color', colors(i,:));
    legend_entries{i} = ['n = ' num2str(n)];
end
plot(theta_max, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
legend(legend_entries, 'Location', 'northeast', 'FontSize', 12);
xlim([0, theta_max+0.2]);
ylim([-0.25, max(f_theta)]);
xticks([0 pi/8 pi/4 3*(pi/8) pi/2]);
xticklabels({'$0$', '$\pi/8$', '$\pi/4$', '$3\pi/8$', '$\pi/2$'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
hold off; grid off;
box on; 
%%

clear;
close all;
clc;
n_values = [0.5, 1, 2]; 
C = 1; 
theta_max = pi/2; 
num_points = 1000;
theta = linspace(0, theta_max, num_points); 

figure;
hold on;
grid on;
xlabel('$\theta$ ($\pi$ radians)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$f_\theta(\theta)$', 'Interpreter', 'latex', 'FontSize', 12);
title('Theta Distribution with Power-Law Decay', 'Interpreter', 'latex', 'FontSize', 14);
legend_entries = cell(length(n_values),1);

colors = lines(length(n_values));
for i = 1:length(n_values)
    n = n_values(i);    
    f_theta = C * ( (pi/2) - theta ).^n;
    f_theta(theta > pi/2) = 0;
    plot(theta, f_theta, 'LineWidth', 2, 'Color', colors(i,:));
    legend_entries{i} = ['n = ' num2str(n)];
end

plot(theta_max, 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
legend(legend_entries, 'Location', 'northeast', 'FontSize', 12);
xlim([0, theta_max]);
ylim([0, max(f_theta)*1.1]);
xticks([0 pi/8 pi/4 3*pi/8 pi/2]);
xticklabels({'$0$', '$\frac{\pi}{8}$', '$\frac{\pi}{4}$', '$\frac{3\pi}{8}$', '$\frac{\pi}{2}$'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
hold off; 
box on;

figure;
hold on;
grid on;
xlabel('$K$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$f_K(K)$', 'Interpreter', 'latex', 'FontSize', 12);
title('Transformed Variable Distribution $K = \frac{1}{\cos(\theta)}$', 'Interpreter', 'latex', 'FontSize', 14);
legend_entries_K = cell(length(n_values),1);


k_max = 10; % Adjust as needed
k = linspace(1, k_max, num_points); 

for i = 1:length(n_values)
    n = n_values(i);    
    f_K = C * k.^-(n+2);
    plot(k, f_K, 'LineWidth', 2, 'Color', colors(i,:));
    legend_entries_K{i} = ['n = ' num2str(n)];
end

legend(legend_entries_K, 'Location', 'northeast', 'FontSize', 12);
xlim([1, k_max]);
ylim([min(C * k.^-(n_values(1)+2)), max(C * k.^-(n_values(end)+2)) * 1.1]);


xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
xticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

hold off; 
box on;
