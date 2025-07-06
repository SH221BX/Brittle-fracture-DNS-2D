%% BASIC of Random picking

clear; clc;

N = 1e6;   R0 = 1.0;     R = 2;     
theta_V = linspace(0, pi/2, N);
f_V = (4/pi) .* (1 - (2*theta_V)/pi);
f_V = f_V / trapz(theta_V, f_V); 
F_V = cumtrapz(theta_V, f_V); 
F_V = F_V / max(F_V);
inverse_CDF = @(u) interp1(F_V, theta_V, u, 'linear', pi/2);
U = rand(N, 1);
theta_V_samples = inverse_CDF(U); 

theta_T_min = 0; theta_T_max = pi/8;
theta_T_samples = theta_T_min + (theta_T_max - theta_T_min) * rand(N, 1);
theta_G_samples = 0 + (pi/2 - 0) * rand(N, 1);

K_V_samples = R0 ./ cos(theta_V_samples); 
K_T_samples = R ./ cos(theta_T_samples);  
K_G_samples = R0 ./ cos(theta_G_samples); 

Z1_samples = min(K_V_samples, K_T_samples);
Z2_samples = min(K_G_samples, K_T_samples); 
X  = 500;
z_min = 1.001; 
z_max = R/cos(theta_T_max);          
z = linspace(z_min, z_max, X);

% f_V = @(theta) (4/pi) .* (1 - (2*theta)/pi); % PDF of V
% f_T = @(theta) (8/pi) .* (theta >= 0 & theta <= pi/8); % PDF of T (uniform)
% f_G = @(theta) (2/pi) .* (theta >= 0 & theta <= pi/2); % PDF of G (uniform)
% 
% K_V = @(theta) R0 ./ cos(theta); % Transformation for K_V
% K_T = @(theta) R ./ cos(theta);  % Transformation for K_T
% K_G = @(theta) R0 ./ cos(theta); % Transformation for K_G
% 
% jacobian = @(k, R) R ./ (k.^2 .* sqrt(1 - (R./k).^2));
% 
% f_KV = @(k) f_V(acos(R0 ./ k)) .* jacobian(k, R0) .* (k >= R0);
% f_KT = @(k) f_T(acos(R ./ k)) .* jacobian(k, R) .* (k >= R & k <= R/cos(pi/8));
% f_KG = @(k) f_G(acos(R0 ./ k)) .* jacobian(k, R0) .* (k >= R0);
% 
% k_min = 1;
% k_max = R / cos(pi/8) - 0.01;
% k = linspace(k_min, k_max, 1000);
% 
% figure('Color', 'w', 'Position', [100, 100, 1200, 400]);
% 
% subplot(1,3,1);
% plot(k, f_KV(k), 'b-', 'LineWidth', 2);
% xlabel('$k$', 'Interpreter', 'latex');
% ylabel('$f_{K_V}(k)$', 'Interpreter', 'latex');
% title('PDF of $K_V$', 'Interpreter', 'latex');
% grid on;
% 
% subplot(1,3,2);
% plot(k, f_KT(k), 'r-', 'LineWidth', 2);
% xlabel('$k$', 'Interpreter', 'latex');
% ylabel('$f_{K_T}(k)$', 'Interpreter', 'latex');
% title('PDF of $K_T$', 'Interpreter', 'latex');
% grid on;
% 
% subplot(1,3,3);
% plot(k, f_KG(k), 'g-', 'LineWidth', 2);
% xlabel('$k$', 'Interpreter', 'latex');
% ylabel('$f_{K_G}(k)$', 'Interpreter', 'latex');
% title('PDF of $K_G$', 'Interpreter', 'latex');
% grid on;


% f_V = @(theta) (4/pi).*(1 - (2.*theta./pi)).*(theta >= 0 & theta <= pi/2);
% f_T = @(theta) (8/pi).*(theta >= 0 & theta <= pi/8);
% f_G = @(theta) (2/pi).*(theta >= 0 & theta <= pi/2);
% 
% f_KV = @(k) f_V(acos(R0./k)) .* (R0./(k.*sqrt(k.^2 - R0.^2))) .* (k >= R0);
% f_KT = @(k) f_T(acos(R./k)) .* (R./(k.*sqrt(k.^2 - R.^2))) .* (k >= R & k <= R./cos(pi/8));
% f_KG = @(k) f_G(acos(R0./k)) .* (R0./(k.*sqrt(k.^2 - R0.^2))) .* (k >= R0);


f_KV = @(z) (4*R0) ./ (pi * z.^2 .* sqrt(1 - (R0./z).^2)) .* (1 - (2/pi)*acos(R0./z)) .* (z >= R0);
f_KT = @(z) (8*R) ./ (pi * z.^2 .* sqrt(1 - (R./z).^2)) .* (z >= R & z <= R/cos(pi/8));
f_KG = @(z) (2*R0) ./ (pi * z.^2 .* sqrt(1 - (R0./z).^2)) .* (z >= R0);

F_KV = @(z) arrayfun(@(zz) integral(@(k) f_KV(k), R0, zz), z);
F_KT = @(z) arrayfun(@(zz) integral(@(k) f_KT(k), R, min(zz, R/cos(pi/8))), z);
F_KG = @(z) arrayfun(@(zz) integral(@(k) f_KG(k), R0, zz), z);

f_Z1 = f_KV(z) .* (1 - F_KT(z)) + f_KT(z) .* (1 - F_KV(z));
f_Z2 = f_KG(z) .* (1 - F_KT(z)) + f_KT(z) .* (1 - F_KG(z));

figure('Color', 'w', 'Position', [100, 100, 1200, 400]);
subplot(1,2,1);
histogram(Z1_samples, 'BinEdges', linspace(z_min, z_max, X), 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 1], 'EdgeColor', 'none');
hold on;
plot(z, f_Z1, '-', 'Color', 'r', 'LineWidth', 2);
xlabel('$z$', 'Interpreter', 'latex');
ylabel('$f_Z(z)$', 'Interpreter', 'latex');
title('$\min(K_V, K_T)$', 'Interpreter', 'latex');
grid on;

subplot(1,2,2);
histogram(Z2_samples, 'BinEdges', linspace(z_min, z_max, X), 'Normalization', 'pdf', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none');
hold on;
plot(z, f_Z2, '-', 'Color', 'b', 'LineWidth', 2);
xlabel('$z$', 'Interpreter', 'latex');
ylabel('$f_Z(z)$', 'Interpreter', 'latex');
title('$\min(K_G, K_T)$', 'Interpreter', 'latex');
grid on;

pV_empirical = mean(K_V_samples < K_T_samples);
pG_empirical = mean(K_G_samples < K_T_samples);
F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(z) f_KT(z), R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, R/cos(pi/8));

lower1 = R0; upper1 = R; lower2 = R; upper2 = R/cos(pi/8);
pG_part1 = integral(f_KG, lower1, upper1);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

fprintf('Empirical P(K_V < K_T) = %.4f\n', pV_empirical);
fprintf('Numerical P(K_V < K_T) = %.4f\n', pV_numerical);
fprintf('Empirical P(K_G < K_T) = %.4f\n', pG_empirical);
fprintf('Numerical P(K_G < K_T) = %.4f\n', pG_numerical);


%%  Random Pick out of two


% R0 = 1 ; R = 1   DONE DONE DONE DONE
clear all; clc; format long;
N = 1e6;
R0 = 1;
R = 1;
n_segments = 30;

pV = zeros(1, n_segments);
P_T = zeros(1, n_segments);
Segments = cell(1, n_segments);
MaxSegments = cell(1, n_segments);

thetaV = linspace(0,pi/2, N);


fVvals = (4/pi) .* (1 - (2 * thetaV ./ pi));
fVvals = fVvals ./ trapz(thetaV, fVvals);
FVvals = cumtrapz(thetaV, fVvals);
FVvals = FVvals ./ max(FVvals);
invFV = @(u) interp1(FVvals, thetaV, u, 'linear', pi/2);

R_max = R/cos(pi/8);
u1 = rand(N,1);
tV1 = invFV(u1);
tT1 = (pi/8) * rand(N,1);
kv1 = R0 ./ cos(tV1);
kv1(kv1 > R_max) = Inf;
kt1 = R ./ cos(tT1);
Z1fresh = min(kv1, kt1);
isV1 = (kv1 < kt1);
pV(1) = mean(isV1);
P_T(1) = 1 - pV(1);
Segments{1} = Z1fresh;
MaxSegments{1} = Segments{1};

for seg = 2:n_segments
    u_seg = rand(N,1);
    tV_seg = invFV(u_seg);
    tT_seg = (pi/8) * rand(N,1);
    kv_seg = R0 ./ cos(tV_seg);
    kv_seg(kv_seg > R_max) = Inf;
    kt_seg = R ./ cos(tT_seg);
    Z1fresh_seg = min(kv_seg, kt_seg);
    isV_seg = (kv_seg < kt_seg);
    tG_seg = (pi/2) * rand(N,1);
    tTb_seg = (pi/8) * rand(N,1);
    kg_seg = R0 ./ cos(tG_seg);
    kv_seg(kg_seg > R_max) = Inf;
    ktb_seg = R ./ cos(tTb_seg);
    Z2fresh_seg = min(kg_seg, ktb_seg);
    isG_seg = (kg_seg < ktb_seg);
    mixFlag = (rand(N,1) < pV(seg-1));
    idxA = randi(N, [N,1]);
    idxB = randi(N, [N,1]);
    Segment = zeros(N,1);
    Segment(mixFlag)  = Z1fresh_seg(idxA(mixFlag));
    Segment(~mixFlag) = Z2fresh_seg(idxB(~mixFlag));
    isT_Seg = false(N,1);
    isT_Seg(mixFlag)  = ~isV_seg(idxA(mixFlag));
    isT_Seg(~mixFlag) = ~isG_seg(idxB(~mixFlag));
    numT = sum(isT_Seg);
    P_T(seg) = numT / N;
    pV(seg) = 1 - P_T(seg);
    Segments{seg} = Segment;
    MaxSegments{seg} = max(MaxSegments{seg-1}, Segments{seg});
end

z_min = 1.00001;
z_max = R_max;
X = 1000;
z = linspace(z_min, z_max, X);

f_KV = @(x) (4 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (1 - (2/pi) * acos(R0 ./ x)) .* (x >= R0 & x <= R / cos(pi/8));
% f_KV = @(x) (2 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));
f_KT = @(x) (8 * R) ./ (pi * x.^2 .* sqrt(1 - (R ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));
f_KG = @(x) (2 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));

F_KV = arrayfun(@(xx) integral(f_KV, R0, xx), z);
F_KT = arrayfun(@(xx) integral(f_KT, R, min(xx, R ./ cos(pi/8))), z);
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
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, R/cos(pi/8));

lower1 = R0; upper1 = R_max; lower2 = R; upper2 = R/cos(pi/8);
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
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

%    f_S{seg} = pV(seg)*fZ1 + (1-pV(seg))*fZ2;
%    F_S{seg} = pV(seg)*F_Z1 + (1-pV(seg))*F_Z2; 

    if seg==1

        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end
end

% figure('Color','w','Position',[100 100 800 400]);

filename = 'R1_two_draw_ran.gif';
fig = figure('Color', 'w');
set(fig, 'Position', [200, 300, 400, 300]);

for seg = 1:n_segments    
    % subplot(5,ceil(n_segments/5),seg);
    histogram(MaxSegments{seg},'BinEdges',linspace(z_min,z_max,50),'Normalization','pdf','FaceColor',[0.7 0.4 0],'EdgeColor','w');
    hold on;
    plot(z,f_Max{seg},'b-','LineWidth',2);
    title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    set(gca, 'FontSize', 14, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

    grid off;
    hold off;
    drawnow; pause(.5);

    % histogram(Segments{seg}, 'BinEdges', linspace(z_min,z_max,X), ...
    %     'Normalization','pdf','FaceColor',[1 0 0],'EdgeColor','none');
    % hold on;
    % plot(z,f_S{seg}, 'k-','LineWidth',1.25);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % grid off;

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if seg == 1 
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else 
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end

end




%%  Random Pick out of two  different ratios


clear all; clc;
format long;
N = 1e6;
R0 = 1; R = 2;
n_segments = 10;

pV = zeros(1,n_segments);
P_T = zeros(1,n_segments);
Segments = cell(1,n_segments);
MaxSegments = cell(1,n_segments);

thetaV = linspace(0,pi/2,N);
thetaV2 = linspace(0,pi/2,N);

fVvals = (4/pi) .* (1 - (2 * thetaV ./ pi));
fVvals = fVvals/trapz(thetaV,fVvals);
FVvals = cumtrapz(thetaV,fVvals);
FVvals = FVvals/max(FVvals);
invFV = @(u) interp1(FVvals,thetaV,u,"spline",pi/2);


u1 = rand(N,1);
tV1 = invFV(u1);
kv1 = R0./cos(tV1);
kv1(kv1 > R) = Inf;
tT1 = (pi/8)*rand(N,1);
kt1 = R./cos(tT1);
Z1fresh = min(kv1,kt1);
isV1 = (kv1 < kt1);
pV(1) = mean(isV1);
P_T(1) = 1 - pV(1);
Segments{1} = Z1fresh;
MaxSegments{1} = Segments{1};


for seg = 2:n_segments
    u_seg = rand(N,1);
    tV_seg = invFV(u_seg);
    kv_seg = R0./cos(tV_seg);
    kv_seg(kv_seg > R) = Inf;
    tT_seg = (pi/8)*rand(N,1);
    kt_seg = R./cos(tT_seg);
    Z1fresh_seg = min(kv_seg,kt_seg);
    isV_seg = (kv_seg < kt_seg);
    tG_seg = (pi/2)*rand(N,1);
    kg_seg = R0./cos(tG_seg);
    kg_seg(kg_seg > R) = Inf;
    tTb_seg = (pi/8)*rand(N,1);
    ktb_seg = R./cos(tTb_seg);
    Z2fresh_seg = min(kg_seg,ktb_seg);
    isG_seg = (kg_seg < ktb_seg);
    mixFlag = (rand(N,1) <= pV(seg-1));
    idxA = randi(N,[N,1]);
    idxB = randi(N,[N,1]);
    Segment = zeros(N,1);
    Segment(mixFlag) = Z1fresh_seg(idxA(mixFlag));
    Segment(~mixFlag) = Z2fresh_seg(idxB(~mixFlag));
    isT_Seg = false(N,1);
    isT_Seg(mixFlag) = ~isV_seg(idxA(mixFlag));
    isT_Seg(~mixFlag) = ~isG_seg(idxB(~mixFlag));
    numT = sum(isT_Seg);
    P_T(seg) = numT/N;
    pV(seg) = 1 - P_T(seg);
    Segments{seg} = Segment;
    MaxSegments{seg} = max(MaxSegments{seg-1},Segment);
    disp(pV(seg));
end

X = 1000;
z_min = 1.0001;
z_max = R/cos(pi/8);
z = linspace(z_min,z_max,X);

f_KV = @(x) (4 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (1 - (2/pi) * acos(R0 ./ x)).*(x>=R0 & x<=R);
% f_KV = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<=R);
f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).*(x>R & x<=R/cos(pi/8));
f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<=R);

F_KV = arrayfun(@(xx) integral(f_KV,R0,xx),z);
F_KT = arrayfun(@(xx) integral(f_KT,R,min(xx,R/cos(pi/8))),z);
F_KG = arrayfun(@(xx) integral(f_KG,R0,xx),z);

F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(z) f_KT(z), R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, R/cos(pi/8));
lower1 = R0; upper1 = R; lower2 = R; upper2 = R/cos(pi/8);
pG_part1 = integral(f_KG, lower1, upper1);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

M = [pV_numerical   pG_numerical
    1-pV_numerical  1-pG_numerical];

P0 = [1; 0];        
           
fZ1 = f_KV(z).*(1-F_KT) + f_KT(z).*(1-F_KV);
fZ2 = f_KG(z).*(1-F_KT) + f_KT(z).*(1-F_KG);
F_Z1 = cumtrapz(z,fZ1); 
F_Z1 = F_Z1/max(F_Z1);
F_Z2 = cumtrapz(z,fZ2); 
F_Z2 = F_Z2/max(F_Z2);

f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
F_Max = cell(1,n_segments);
f_Max = cell(1,n_segments);
disp('______________________');

for seg = 1:n_segments
    if seg == 1
        P_N = P0;
    else
        [V, D] = eig(M);
        D(1,1) = D(1,1)^(seg-1);
        D(2,2) = D(2,2)^(seg-1);
        M_pow = V * D / V;
        P_N = M_pow * P0;
        disp(P_N(1)); 
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
end

% figure('Color','w','Position',[100 100 800 400]);

% filename = 'R12_two_draw_ran.gif';
% fig = figure('Color', 'w','MenuBar','none');
% set(fig, 'Position', [200, 300, 400, 300]);

for seg = 1:n_segments    
    % subplot(5,ceil(n_segments/5),seg);
    % histogram(MaxSegments{seg},'BinEdges',linspace(z_min,z_max,50),'Normalization','pdf','FaceColor',[0.7 0.4 0],'EdgeColor','w');
    % hold on;
    % plot(z,f_Max{seg},'b-','LineWidth',1.5);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 14, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % 
    % grid off;
    % hold off;
    % drawnow; pause(0);

    % histogram(Segments{seg}, 'BinEdges', linspace(z_min,z_max,X), ...
    %     'Normalization','pdf','FaceColor',[1 0 0],'EdgeColor','none');
    % hold on;
    % plot(z,f_S{seg}, 'k-','LineWidth',1.25);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % grid off;

    % frame = getframe(gcf);
    % img = frame2im(frame);
    % [imind, cm] = rgb2ind(img, 256);
    % 
    % if seg == 1 
    %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    % else 
    %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    % end

    norm_factor = trapz(z,f_Max{seg});
    mean_pdf(seg) = trapz(z,z.*f_Max{seg})/norm_factor;
end

figure;
plot(1:n_segments,mean_pdf,'c-','LineWidth',2);
xlabel('N','Interpreter','latex','FontSize',14);
ylabel('$\mathrm{E[K_{M}]}$','Interpreter','latex','FontSize',14);
grid on;



%% VORONOI

% R0 = 1 ; R = 1   DONE DONE DONE DONE
clear; clc; format long;
N = 1e6;
R0 = 1;
R = 1;
n_segments = 30;

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

R_max = R/cos(pi/8);
u1 = rand(N,1);
tV1 = invFV(u1);
tT1 = (pi/8) * rand(N,1);
kv1 = R0 ./ cos(tV1);
kv1(kv1 > R_max) = Inf;
kt1 = R ./ cos(tT1);
Z1fresh = min(kv1, kt1);
isV1 = (kv1 < kt1);
pV(1) = mean(isV1);
P_T(1) = 1 - pV(1);
Segments{1} = Z1fresh;
MaxSegments{1} = Segments{1};

for seg = 2:n_segments
    u_seg = rand(N,1);
    tV_seg = invFV(u_seg);
    tT_seg = (pi/8) * rand(N,1);
    kv_seg = R0 ./ cos(tV_seg);
    kv_seg(kv_seg > R_max) = Inf;
    kt_seg = R ./ cos(tT_seg);
    Z1fresh_seg = min(kv_seg, kt_seg);
    isV_seg = (kv_seg < kt_seg);
    tG_seg = (pi/2) * rand(N,1);
    tTb_seg = (pi/8) * rand(N,1);
    kg_seg = R0 ./ cos(tG_seg);
    kv_seg(kg_seg > R_max) = Inf;
    ktb_seg = R ./ cos(tTb_seg);
    Z2fresh_seg = min(kg_seg, ktb_seg);
    isG_seg = (kg_seg < ktb_seg);
    mixFlag = (rand(N,1) < pV(seg-1));
    idxA = randi(N, [N,1]);
    idxB = randi(N, [N,1]);
    Segment = zeros(N,1);
    Segment(mixFlag)  = Z1fresh_seg(idxA(mixFlag));
    Segment(~mixFlag) = Z2fresh_seg(idxB(~mixFlag));
    isT_Seg = false(N,1);
    isT_Seg(mixFlag)  = ~isV_seg(idxA(mixFlag));
    isT_Seg(~mixFlag) = ~isG_seg(idxB(~mixFlag));
    numT = sum(isT_Seg);
    P_T(seg) = numT / N;
    pV(seg) = 1 - P_T(seg);
    Segments{seg} = Segment;
    MaxSegments{seg} = max(MaxSegments{seg-1}, Segments{seg});
end

z_min = 1.00001;
z_max = R_max;
X = 1000;
z = linspace(z_min, z_max, X);

f_KV = @(x) (A ./ (1 + exp(M * ((pi/2 - acos(R0 ./x)).^(1/3) - phi))) ) .* (R0 ./ x.^2) ./ sqrt(1 - (R0 ./x).^2) .* (x >= R & x <= R / cos(pi/8));
f_KT = @(x) (8 * R) ./ (pi * x.^2 .* sqrt(1 - (R ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));
f_KG = @(x) (2 * R0) ./ (pi * x.^2 .* sqrt(1 - (R0 ./ x).^2)) .* (x >= R & x <= R / cos(pi/8));

F_KV = arrayfun(@(xx) integral(f_KV, R0, xx), z);
F_KT = arrayfun(@(xx) integral(f_KT, R, min(xx, R ./ cos(pi/8))), z);
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
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, R/cos(pi/8));

lower1 = R0; upper1 = R_max; lower2 = R; upper2 = R/cos(pi/8);
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
        disp(P_N);
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

%    f_S{seg} = pV(seg)*fZ1 + (1-pV(seg))*fZ2;
%    F_S{seg} = pV(seg)*F_Z1 + (1-pV(seg))*F_Z2; 

    if seg==1
        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end
end

filename = 'R1_two_draw_vor.gif';
fig = figure('Color', 'w','MenuBar','none');
set(fig, 'Position', [200, 300, 400, 300]);

for seg = 1:n_segments    
    % subplot(5,ceil(n_segments/5),seg);
    histogram(MaxSegments{seg},'BinEdges',linspace(z_min,z_max,50),'Normalization','pdf','FaceColor',[0.4 0.8 0.5],'EdgeColor','w');
    hold on;
    plot(z,f_Max{seg},'r-','LineWidth',1.5);
    title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    set(gca, 'FontSize', 14, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');

    grid off;
    hold off;
    drawnow; pause(0);

    % histogram(Segments{seg}, 'BinEdges', linspace(z_min,z_max,X), ...
    %     'Normalization','pdf','FaceColor',[1 0 0],'EdgeColor','none');
    % hold on;
    % plot(z,f_S{seg}, 'k-','LineWidth',1.25);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % grid off;

    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);

    if seg == 1 
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    else 
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end

end


%%   VORONOI with KGB/KT ratio 


clear all; clc; format long;
N = 1e6;
R0 = 1; R = 2;
n_segments = 100;

pV = zeros(1,n_segments);
P_T = zeros(1,n_segments);
Segments = cell(1,n_segments);
MaxSegments = cell(1,n_segments);

thetaV = linspace(0,pi/2,N);
thetaV2 = linspace(0,pi/2,N);

% A = 1.09530; M = -12.2242; phi = 0.8483;
A = 1.052;  M = -12.09;  phi = 0.825;
fV_unnormalized = 1./(1+exp(M*((pi/2 - thetaV).^(1/3)-phi)));
C = trapz(thetaV,fV_unnormalized);
A = 1/C;


fVvals = A./(1+exp(M*((pi/2 - thetaV).^(1/3)-phi)));
fVvals = fVvals/trapz(thetaV,fVvals);
FVvals = cumtrapz(thetaV,fVvals);
FVvals = FVvals/max(FVvals);
invFV = @(u) interp1(FVvals,thetaV,u,"spline",pi/2);


u1 = rand(N,1);
tV1 = invFV(u1);
kv1 = R0./cos(tV1);
kv1(kv1 > R) = Inf;
tT1 = (pi/8)*rand(N,1);
kt1 = R./cos(tT1);
Z1fresh = min(kv1,kt1);
isV1 = (kv1 < kt1);
pV(1) = mean(isV1);
P_T(1) = 1 - pV(1);
Segments{1} = Z1fresh;
MaxSegments{1} = Segments{1};


for seg = 2:n_segments
    u_seg = rand(N,1);
    tV_seg = invFV(u_seg);
    kv_seg = R0./cos(tV_seg);
    kv_seg(kv_seg > R) = Inf;
    tT_seg = (pi/8)*rand(N,1);
    kt_seg = R./cos(tT_seg);
    Z1fresh_seg = min(kv_seg,kt_seg);
    isV_seg = (kv_seg < kt_seg);
    tG_seg = (pi/2)*rand(N,1);
    kg_seg = R0./cos(tG_seg);
    kg_seg(kg_seg > R) = Inf;
    tTb_seg = (pi/8)*rand(N,1);
    ktb_seg = R./cos(tTb_seg);
    Z2fresh_seg = min(kg_seg,ktb_seg);
    isG_seg = (kg_seg < ktb_seg);
    mixFlag = (rand(N,1) <= pV(seg-1));
    idxA = randi(N,[N,1]);
    idxB = randi(N,[N,1]);
    Segment = zeros(N,1);
    Segment(mixFlag) = Z1fresh_seg(idxA(mixFlag));
    Segment(~mixFlag) = Z2fresh_seg(idxB(~mixFlag));
    isT_Seg = false(N,1);
    isT_Seg(mixFlag) = ~isV_seg(idxA(mixFlag));
    isT_Seg(~mixFlag) = ~isG_seg(idxB(~mixFlag));
    numT = sum(isT_Seg);
    P_T(seg) = numT/N;
    pV(seg) = 1 - P_T(seg);
    Segments{seg} = Segment;
    MaxSegments{seg} = max(MaxSegments{seg-1},Segment);
%   disp(pV(seg));
end

X = 1000;
z_min = 1.0001;
z_max = R/cos(pi/8);
z = linspace(z_min,z_max,X);

f_KV = @(x) (A./(1+exp(M*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).*(x>=R0 & x<=R);
f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).*(x>R & x<=R/cos(pi/8));
f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<=R);

F_KV = arrayfun(@(xx) integral(f_KV,R0,xx),z);
F_KT = arrayfun(@(xx) integral(f_KT,R,min(xx,R/cos(pi/8))),z);
F_KG = arrayfun(@(xx) integral(f_KG,R0,xx),z);

F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(z) f_KT(z), R, kk), k);
pV_numerical = integral(@(k) f_KV(k) .* F_KT_survival(k), R0, R/cos(pi/8));

lower1 = R0; upper1 = R; lower2 = R; upper2 = R/cos(pi/8);
pG_part1 = integral(f_KG, lower1, upper1);
pG_part2 = integral(@(k) f_KG(k) .* F_KT_survival(k), lower2, upper2);
pG_numerical = pG_part1 + pG_part2;

M = [pV_numerical   pG_numerical
    1-pV_numerical  1-pG_numerical];
P0 = [1; 0];        
           
fZ1 = f_KV(z).*(1-F_KT) + f_KT(z).*(1-F_KV);
fZ2 = f_KG(z).*(1-F_KT) + f_KT(z).*(1-F_KG);
F_Z1 = cumtrapz(z,fZ1); 
F_Z1 = F_Z1/max(F_Z1);
F_Z2 = cumtrapz(z,fZ2); 
F_Z2 = F_Z2/max(F_Z2);

f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
F_Max = cell(1,n_segments);
f_Max = cell(1,n_segments);


for seg = 1:n_segments
    if seg == 1
        P_N = P0;
    else
        [V, D] = eig(M);
        D(1,1) = D(1,1)^(seg-1);
        D(2,2) = D(2,2)^(seg-1);
        M_pow = V * D / V;
        P_N = M_pow * P0;
        disp(P_N);
    end

    f_S{seg} = P_N(1)*fZ1 + (1-P_N(1))*fZ2;
    F_S{seg} = P_N(1)*F_Z1 + (1-P_N(1))*F_Z2;

%   f_S{seg} = pV(seg)*fZ1 + (1-pV(seg))*fZ2;
%   F_S{seg} = pV(seg)*F_Z1 + (1-pV(seg))*F_Z2; 

    if seg==1
        F_Max{seg} = F_S{seg};
        f_Max{seg} = f_S{seg};
    else
        F_Max{seg} = F_Max{seg-1}.*F_S{seg};
        f_Max{seg} = f_Max{seg-1}.*F_S{seg} + F_Max{seg-1}.*f_S{seg};
    end
end

% filename = 'R12_vor.gif';
% fig = figure('Color', 'w','MenuBar','none');
% set(fig, 'Position', [200, 300, 400, 300]);

for seg = 1:n_segments    
    % subplot(5,ceil(n_segments/5),seg);
    % histogram(MaxSegments{seg},'BinEdges',linspace(z_min,z_max,50),'Normalization','pdf','FaceColor',[0.4 0.8 0.5],'EdgeColor','w');
    % hold on;
    % plot(z,f_Max{seg},'r-','LineWidth',1.5);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 14, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % 
    % grid off; hold off;
    % drawnow; pause(0);

    % histogram(Segments{seg}, 'BinEdges', linspace(z_min,z_max,X), ...
    %     'Normalization','pdf','FaceColor',[1 0 0],'EdgeColor','none');
    % hold on;
    % plot(z,f_S{seg}, 'k-','LineWidth',1.25);
    % title(['$\rm N = \, ' num2str(seg) '$'],'Interpreter','latex');
    % set(gca, 'FontSize', 8, 'LineWidth', 1, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
    % grid off;

    % frame = getframe(gcf);
    % img = frame2im(frame);
    % [imind, cm] = rgb2ind(img, 256);

    % if seg == 1 
    %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
    % else 
    %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    % end
    norm_factor = trapz(z,f_Max{seg});
    mean_pdf(seg) = trapz(z,z.*f_Max{seg})/norm_factor;
end

% figure;
plot(1:n_segments,mean_pdf,'b-','LineWidth',2);
% xlabel('N','Interpreter','latex','FontSize',14);
% ylabel('$\mathrm{E[K_{M}]}$','Interpreter','latex','FontSize',14);
% grid on;

%%


clear all;


A = 1.0530; mval = -12.08; phi = 0.8483;
R0 = 1;

for i= 1.001:0.1:3
R = i;
n_segments = 40;
X = 1000;
z_min = 1.0001;
z_max = R/cos(pi/8);
z = linspace(z_min,z_max,X);
f_KV = @(x) (A./(1+exp(mval*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).*(x>=R0 & x<=R);
f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).*(x>R & x<=R/cos(pi/8));
f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<=R);
F_KV = arrayfun(@(xx) integral(f_KV,R0,xx),z);
F_KT = arrayfun(@(xx) integral(f_KT,R,min(xx,R/cos(pi/8))),z);
F_KG = arrayfun(@(xx) integral(f_KG,R0,xx),z);
F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(zz) f_KT(zz),R,kk),k);
pV_numerical = integral(@(k) f_KV(k).*F_KT_survival(k),R0,R/cos(pi/8));
lower1 = R0;
upper1 = R;
lower2 = R;
upper2 = R/cos(pi/8);
pG_part1 = integral(f_KG,lower1,upper1);
pG_part2 = integral(@(k) f_KG(k).*F_KT_survival(k),lower2,upper2);
pG_numerical = pG_part1 + pG_part2;
M = [pV_numerical pG_numerical; 1-pV_numerical 1-pG_numerical];
P0 = [1;0];
fZ1 = f_KV(z).*(1-F_KT) + f_KT(z).*(1-F_KV);
fZ2 = f_KG(z).*(1-F_KT) + f_KT(z).*(1-F_KG);
F_Z1 = cumtrapz(z,fZ1);
F_Z1 = F_Z1/max(F_Z1);
F_Z2 = cumtrapz(z,fZ2);
F_Z2 = F_Z2/max(F_Z2);
pT_V = trapz(z,f_KT(z).*(1-F_KV))./trapz(z,fZ1);
pT_G = trapz(z,f_KT(z).*(1-F_KG))./trapz(z,fZ2);
f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
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
    T_count = seg*(P_N(1)*pT_V + P_N(2)*pT_G);

if seg == 40
    loglog(i,(P_N(1)*pT_V + P_N(2)*pT_G),'ro'); hold on;
end
     % plot(seg,T_count,'ro'); hold on;
     
end

end

%%

clear all; 
A = 1.0530; mval = -12.08; phi = 0.83;

 
R0 = 1;
R = 2;
n_segments = 30;
X = 1000;
z_min = 1.0001;
z_max = R/cos(pi/8);
z = linspace(z_min,z_max,X);
f_KV = @(x) (A./(1+exp(mval*((pi/2 - acos(R0./x)).^(1/3)-phi)))).*(R0./x.^2)./sqrt(1-(R0./x).^2).*(x>=R0 & x<=R);
f_KT = @(x) (8*R)./(pi*x.^2.*sqrt(1-(R./x).^2)).*(x>R & x<=R/cos(pi/8));
f_KG = @(x) (2*R0)./(pi*x.^2.*sqrt(1-(R0./x).^2)).*(x>=R0 & x<=R);
F_KV = arrayfun(@(xx) integral(f_KV,R0,xx),z);
F_KT = arrayfun(@(xx) integral(f_KT,R,min(xx,R/cos(pi/8))),z);
F_KG = arrayfun(@(xx) integral(f_KG,R0,xx),z);
F_KT_survival = @(k) 1 - arrayfun(@(kk) integral(@(zz) f_KT(zz),R,kk),k);
pV_numerical = integral(@(k) f_KV(k).*F_KT_survival(k),R0,R/cos(pi/8));
lower1 = R0;
upper1 = R;
lower2 = R;
upper2 = R/cos(pi/8);
pG_part1 = integral(f_KG,lower1,upper1);
pG_part2 = integral(@(k) f_KG(k).*F_KT_survival(k),lower2,upper2);
pG_numerical = pG_part1 + pG_part2;
M = [pV_numerical pG_numerical; 1-pV_numerical 1-pG_numerical];
P0 = [1;0];
fZ1 = f_KV(z).*(1-F_KT) + f_KT(z).*(1-F_KV);
fZ2 = f_KG(z).*(1-F_KT) + f_KT(z).*(1-F_KG);
F_Z1 = cumtrapz(z,fZ1);
F_Z1 = F_Z1/max(F_Z1);
F_Z2 = cumtrapz(z,fZ2);
F_Z2 = F_Z2/max(F_Z2);

pT_V = trapz(z,f_KT(z).*(1-F_KV))./trapz(z,fZ1);
pT_G = trapz(z,f_KT(z).*(1-F_KG))./trapz(z,fZ2);

f_S = cell(1,n_segments);
F_S = cell(1,n_segments);
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
    T_count = seg*(P_N(1)*pT_V + P_N(2)*pT_G);
    disp(T_count);
end

% pt = [pT_V; pT_G];
% dist3d = zeros(n_segments+1,n_segments+1,2);
% dist3d(1,1,1) = 1;
% for seg = 1:n_segments
%     for k = 0:(seg-1)
%         for st = 1:2
%             p = dist3d(seg,k+1,st);
%             if p>0
%                 for j = 1:2
%                     dist3d(seg+1,k+1,j) = dist3d(seg+1,k+1,j) + p*(1-pt(st))*M(st,j);
%                     dist3d(seg+1,k+2,j) = dist3d(seg+1,k+2,j) + p*pt(st)*M(st,j);
%                 end
%             end
%         end
%     end
% end
% figure;
% rows = 5; cols = 10;
% for seg = 1:n_segments
%     if seg == 30
%     v = sum(dist3d(seg+1,1:seg+1,:),3);
%     % subplot(rows,cols,seg);
%     bar(0:seg,v(1:seg+1));
%     xlim([-0.5 seg+0.5]);
%     ylim([0 1]);
%     end
% end







