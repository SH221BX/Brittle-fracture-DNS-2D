clear; clc; close all;
E_AA = 1;
a = 1;
beta = linspace(0, pi/2, 200);
gamma_beta = (E_AA/(1*a))*(cos(beta) + sin(beta));
figure; plot(beta, gamma_beta, 'LineWidth',2);
xlabel('\beta','Interpreter','tex');
ylabel('\gamma(\beta)','Interpreter','tex');
% title('Surface Energy vs. \theta','Interpreter','tex');
% grid on;

%%
clear; clc; close all;
E_AA = 1; a = 1;
beta = linspace(0, pi/2, 200);
gamma_beta = (E_AA/(a))*(cos(beta) + sin(beta));

figure;
polarplot(beta, gamma_beta, 'LineWidth',2);
title('Surface Energy vs. \beta','Interpreter','tex');

%%  cusp 90 deg apart

clear; clc; close all;

E_AA = 1;  a    = 1;         
theta = linspace(0, 2*pi, 400);
gamma_val = (E_AA/a) * ( abs(cos(theta)) + abs(sin(theta)) );
figure;
polarplot(theta, gamma_val, 'LineWidth',1);

%% Multiple even cusp 45 deg apart

clear; clc; close all;

E_AA = 1;   a    = 1;         
theta = linspace(0, 2*pi, 400);
gamma_val = (E_AA/a) * (abs(cos(4*theta)) + abs(sin(4*theta)) );

figure;
polarplot(theta, gamma_val, 'LineWidth',2);

%%
clear; clc; close all;

E_AA = 1;  
a    = 1;         
theta = linspace(0, 2*pi, 400);

w_cos = 1;    
w_sin = 1.25;  

gamma_val = (E_AA/a) * ( w_cos * abs(cos(2*theta)) + w_sin * abs(sin(2*theta)) );

figure;
polarplot(theta, gamma_val, 'LineWidth', 2);



%% gamma constant (isotropic surface energy)

clear; clc; close all;

E_AA = 1;  a    = 1;         
theta = linspace(0, 2*pi, 400);
gamma_val = (E_AA/a) * ones(size(theta));
figure;
polarplot(theta, gamma_val, 'LineWidth',2);

%%
clear; clc; close all;
gamma0 = 1;        
M = 1e6;          
epsilon = 1e-2;    
theta = linspace(0, 2*pi, 400);

gamma_val = M * ones(size(theta));

idx = ( abs(theta) < epsilon | ...
        abs(theta - 2*pi) < epsilon | ...
        abs(theta - pi/2) < epsilon | ...
        abs(theta - pi) < epsilon | ...
        abs(theta - 3*pi/2) < epsilon );

gamma_val(idx) = gamma0;
figure;
polarplot(theta, gamma_val, 'LineWidth', 1);
title('Surface Energy with Cusps at [1,0] and [0,1] Directions');

%%
clear; clc; close all;

gamma0 = 1;       
M = 1.5;          
epsilon = 1e-2;    

theta = linspace(0, 2*pi, 1000);
gamma_val = M * ones(size(theta));

idx = ( abs(theta) < epsilon | ...
        abs(theta - 2*pi) < epsilon | ...
        abs(theta - pi/2) < epsilon | ...
        abs(theta - pi) < epsilon | ...
        abs(theta - 3*pi/2) < epsilon );

gamma_val(idx) = gamma0;
figure;
polarplot(theta, gamma_val,'Color','r' ,'LineWidth',2);




%%

clear; clc; close all;
E_AA = 1;
a = 1;
theta = linspace(0, 2*pi, 400);
gamma_val = (E_AA/a) * (abs(cos(theta)) + abs(sin(theta)));

figure;
for phi = linspace(0,1,1)
    rotated_theta = theta + phi;
    polarplot(rotated_theta, gamma_val, 'LineWidth', 2);
    title(['Rotated Gamma System by ', num2str(phi*180/pi, '%.1f'), '°'], 'Interpreter','tex');
    drawnow;
    pause(0.1);
end


%%

clear; clc; close all;
E_AA = 1;
a = 1;
theta = linspace(0, 2*pi, 400);
gamma_val = (E_AA/a) * ( abs(cos(theta)) + abs(sin(theta)) );

figure; 
pax = polaraxes;
hold(pax, 'on');
r_max = max(gamma_val); 

nSteps = 1; % Number of rotation steps

for phi = linspace(0, pi/6, nSteps)

    polarplot(pax, theta + phi, gamma_val, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    r_line = linspace(-r_max, r_max, 100);
    theta_line = phi * ones(size(r_line));
    polarplot(pax, theta_line, r_line, 'r--', 'LineWidth', 3);
    polarplot(pax, phi+pi/2 * ones(size(r_line)), r_line, 'r--', 'LineWidth', 3);


    polarplot(pax,  0 * ones(size(r_line)), r_line, 'k-', 'LineWidth', 2);
     polarplot(pax,  pi/2 * ones(size(r_line)), r_line, 'k-', 'LineWidth', 2);
    
    % title(pax, ['Rotated Gamma System by ', num2str(phi*180/pi, '%.1f'), '°'], 'Interpreter', 'tex');
    drawnow;
    % pause(0.1);
    % cla(pax); % Clear the axes for the next frame
    % hold(pax, 'on');
end


%%

clear; clc; close all;

E_AA = 1; a = 1;
theta = linspace(0, pi/2, 200);
gamma_theta = (E_AA / (1*a)) * (abs(cos(theta)) + abs(sin(theta)));

theta_g = pi/6;  
gamma_theta_rot = (E_AA / a) * (abs(cos(theta + theta_g)) + abs(sin(theta + theta_g)));

figure;
plot(theta, gamma_theta, 'LineWidth', 2); hold on;
plot(theta, gamma_theta_rot, 'LineWidth', 2); 
xlabel('\beta', 'Interpreter', 'tex');
ylabel('\gamma(\beta)', 'Interpreter', 'tex');

legend('$\theta = 0$','$\theta = \frac{\pi}{6}$','Location','Best');
grid on;

% figure;
% polarplot(theta, gamma_theta_rot, 'LineWidth',1);

%%

clear; clc; close all;

E_AA   = 1;
a      = 1;
theta  = linspace(0, pi/4, 200); 
theta_g = pi/4;                   

gamma_theta = (E_AA/a) * (cos(theta) + sin(theta));
gamma_theta_rot = (E_AA/a) * (cos(theta + theta_g) + sin(theta + theta_g));

K_original = sqrt(gamma_theta)      ./ cos(theta);
K_rotated  = sqrt(gamma_theta_rot)  ./ cos(theta);

figure;
plot(theta, K_original, 'LineWidth',2); hold on;
plot(theta, K_rotated,  'LineWidth',2);
xlabel('\theta', 'Interpreter', 'tex');
ylabel('K(\theta)', 'Interpreter', 'tex');
% legend('K_{original}(\theta)','K_{rotated}(\theta + \theta_g)','Location','best');
grid on;
legend('$\theta = 0$','$\theta = \frac{\pi}{6}$','Location','Best');

%%

clear; clc; close all;
E_AA = 1;
a    = 1;
theta = linspace(0, pi/4, 200);                               % theta values along x-axis
theta_g_list = [0 pi/6 pi/4];    %linspace(0, pi/4-0.2, 5);        % 5 equally spaced theta_g values from 0 to pi/4

figure; hold on;
legend_entries = cell(length(theta_g_list), 1);


for i = 1:length(theta_g_list)
    theta_g = theta_g_list(i);
    
    gamma_theta_rot = (E_AA/a) * (abs(cos(theta + theta_g)) + abs(sin(theta + theta_g)));
    K = sqrt(gamma_theta_rot) ./ abs(cos(theta + theta_g));
   
    plot(theta, K, 'LineWidth', 2);
    legend_entries{i} = sprintf('\\theta_g = %.2f', theta_g);
end

ylim([1 5]);
xlabel('\theta', 'Interpreter', 'tex');
ylabel('K', 'Interpreter', 'tex');
legend(legend_entries, 'Location', 'best');
legend('$\theta = 0$','$\theta = \frac{\pi}{6}$','Location','Best');
grid on;


%%
clear all;
N = 1e6;   
rng('default');

Gamma = 1 + rand(N,1);
Theta = (pi/2)*rand(N,1);

figure;

subplot(1,2,1);
histogram(Gamma,'Normalization','pdf');
title('Distribution of \Gamma \sim Uniform(1,2)');
xlabel('\Gamma');
ylabel('Density');

subplot(1,2,2);
histogram(Theta,'Normalization','pdf');
title('Distribution of \Theta \sim Uniform(0,\pi/2)');
xlabel('\Theta');
ylabel('Density');


%%

clear all;
N = 1e7;
rng('default');
g = 1 + rand(N,1);
t = (pi/2)*rand(N,1);
K = sqrt(g)./cos(t);

kmin = 1;
kmax = prctile(K,99.5);
edges = linspace(kmin,kmax,10000);
[count,edges] = histcounts(K,edges,'Normalization','pdf');
centers = 0.5*(edges(1:end-1)+edges(2:end));

pdf_vals = arrayfun(@(x) pdfK(x), centers);

bar(centers,count,'hist')
hold on;
plot(centers,pdf_vals,'r-','LineWidth',2)
hold off;
xlim([1 10]);
xlabel('K', 'Interpreter', 'latex');
ylabel('$f_{K}$', 'Interpreter', 'latex');
%%

clear all;
N = 1e7;
Q = 100;  
rng('default');
g = 1 + (Q-1)*rand(N,1);
t = (pi/2)*rand(N,1);
K = sqrt(g)./cos(t);

kmin = 1;
kmax = prctile(K,99.5);
edges = linspace(kmin,kmax,10000);
[count,edges] = histcounts(K,edges,'Normalization','pdf');
centers = 0.5*(edges(1:end-1)+edges(2:end));

pdf_vals = zeros(size(centers));
idx1 = centers < sqrt(Q);      % Case: 1 <= k < sqrt(Q)
idx2 = centers >= sqrt(Q);     % Case: k >= sqrt(Q)

theta_a = acos(1 ./ centers(idx1));
F_a = theta_a/2 + sin(2*theta_a)/4;
pdf_vals(idx1) = (4 .* centers(idx1) ./ (pi*(Q-1))) .* F_a;

theta1 = acos(1 ./ centers(idx2));
theta2 = acos(sqrt(Q) ./ centers(idx2));
F1 = theta1/2 + sin(2*theta1)/4;
F2 = theta2/2 + sin(2*theta2)/4;
pdf_vals(idx2) = (4 .* centers(idx2) ./ (pi*(Q-1))) .* (F1 - F2);

bar(centers, count, 'hist');
hold on;
plot(centers, pdf_vals, 'r-', 'LineWidth',2);
hold off;
xlim([1 50]);
xlabel('K','Interpreter','latex');
ylabel('$f_{K}$','Interpreter','latex');
title(['$Q = ', num2str(Q), '$'],'Interpreter','latex');


%%
function val = pdfK(x)
if x < 1
    val = 0;
elseif x < sqrt(2)
    val = (2/pi)*( x*acos(1/x) + sqrt(x^2 - 1)/x );
else
    val = (2/pi)*( x*( acos(1/x) - acos(sqrt(2)/x) ) + ...
                   ( sqrt(x^2 - 1) - sqrt(2)*sqrt(x^2 - 2 ) )/x );
end
end




