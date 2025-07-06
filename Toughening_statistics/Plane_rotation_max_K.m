
% Plane 4

clc;clear all;

theta = -22.5:0.1:22.5;
p = 1;

y1 = abs(1 ./ cosd(theta).^p);             % Plane 1
y2 = abs(1 ./ cosd(theta + 90).^p);        % Plane 3
y3 = abs(1 ./ cosd(theta + 45).^p);        % Plane 2
y4 = abs(1 ./ cosd(theta + 135).^p);       % Plane 4

combined_y = [y1; y2; y3; y4];
[min_values, min_indices] = min(combined_y, [], 1); 
theta_min = theta(min_indices);

figure;
hold on;
plot(theta, y1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Plane 1: $1/\cos^2(\theta)$');
plot(theta, y2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Plane 3: $1/\cos^2(\theta + 90)$');
plot(theta, y3, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Plane 2: $1/\cos^2(\theta + 45)$');
plot(theta, y4, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Plane 4: $1/\cos^2(\theta - 45)$');

grid off;
ylabel('K_{Ic}', 'Interpreter', 'tex');
xlabel('Theta (degrees)', 'Interpreter', 'latex');
ylim([1 5]);box on;

 figure;
 plot(theta, min_values, 'r-', 'MarkerSize', 5, 'DisplayName', 'Min of all planes','LineWidth', 1.5); hold on;
 xlabel('Theta (degrees)', 'Interpreter', 'latex');
 ylabel('K_{Ic} Minimum', 'Interpreter', 'tex');

grid off;
xlim([-90 90]);
ylim([1 1.2]);  
 box on;
disp(max(min_values));
%% Plane 03


clc;clear all;

theta = -30:0.1:30;
p = 1;

y1 = abs(1 ./ cosd(theta).^p);             % Plane 1
y2 = abs(1 ./ cosd(theta + 60).^p);        % Plane 3
y3 = abs(1 ./ cosd(theta - 60).^p);       % Plane 2


combined_y = [y1; y2; y3];
[min_values, min_indices] = min(combined_y, [], 1); 
theta_min = theta(min_indices);

figure;
hold on;
plot(theta, y1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Plane 1: $1/\cos^2(\theta)$');
plot(theta, y2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Plane 3: $1/\cos^2(\theta + 90)$');
plot(theta, y3, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Plane 2: $1/\cos^2(\theta + 45)$');


grid off;
ylabel('K_{Ic}', 'Interpreter', 'tex');
xlabel('Theta (degrees)', 'Interpreter', 'latex');
ylim([1 5]);box on;

 figure;
 plot(theta, min_values, 'b-', 'MarkerSize', 5, 'DisplayName', 'Min of all planes','LineWidth', 1.5); hold on;
 xlabel('Theta (degrees)', 'Interpreter', 'latex');
 ylabel('K_{Ic} Minimum', 'Interpreter', 'tex');

grid off;
xlim([-90 90]);
ylim([1 1.2]);  
 box on;
disp(max(min_values));


%% Plane_02


clc;clear all;

theta = -45:0.5:45;
p = 1;

y1 = abs(1 ./ cosd(theta).^p);            % Plane 1
y2 = abs(1 ./ cosd(theta + 90).^p);       % Plane 2


combined_y = [y1; y2];
[min_values, min_indices] = min(combined_y, [], 1); 
theta_min = theta(min_indices);

figure;
hold on;
plot(theta, y1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Plane 1: $1/\cos^2(\theta)$');
plot(theta, y2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Plane 3: $1/\cos^2(\theta + 90)$');

grid off;
ylabel('K_{Ic}', 'Interpreter', 'tex');
xlabel('Theta (degrees)', 'Interpreter', 'latex');
ylim([1 5]);box on;

 figure;
 plot(theta, min_values, 'c-', 'MarkerSize', 5, 'DisplayName', 'Min of all planes','LineWidth', 1.5);
 xlabel('Theta (degrees)', 'Interpreter', 'latex');
 ylabel('K_{Ic} Minimum', 'Interpreter', 'tex');
disp(max(min_values));
% grid off;
% xlim([-90 90]);
% ylim([1 1.2]);  
%  box on;

%% Plane_01

clc;clear all;
theta = -90:0.5:90;
p = 1;

y1 = abs(1 ./ cosd(theta).^p);            % Plane 1


combined_y = [y1];
[min_values, min_indices] = min(combined_y, [], 1); 
theta_min = theta(min_indices);

figure;
hold on;
plot(theta, y1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Plane 1: $1/\cos^2(\theta)$');

grid off;
ylabel('K_{Ic}', 'Interpreter', 'tex');
xlabel('Theta (degrees)', 'Interpreter', 'latex');
ylim([1 5]);box on;

 figure;
 plot(theta, min_values, 'r-', 'MarkerSize', 5, 'DisplayName', 'Min of all planes','LineWidth', 1.5);
 xlabel('Theta (degrees)', 'Interpreter', 'latex');
 ylabel('K_{Ic} Minimum', 'Interpreter', 'tex');














