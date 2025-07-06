clear all; clc
N = 1e7;  K0 = 1;
planes = [];
planes1 = [1  0  0;
          -1  0  0;
           0  1  0;
           0 -1  0;
           0  0  1;
           0  0 -1]';
planes = [planes  planes1];

minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C = abs(R(3,:)*planes);
    K_eq = K0./max(eps,C);
    minKeq_vals(i) = min(K_eq);
end

figure;
subplot(1,2,1);
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','c','FaceAlpha',1,'EdgeColor','c');hold on;


theta1 = linspace(0, pi/4, 1000);
theta2 = linspace(pi/4, atan(sqrt(2)), 100);

pdf1 = 3 * sin(theta1);
pdf2 = 3 *(2/pi)* (asin(cot(theta2)) - acos(cot(theta2))) .* sin(theta2);
plot(theta1, pdf1, 'r', theta2, pdf2, 'r','LineWidth',2)
ylim([0 2.5])
xlabel('$\theta$','Interpreter','latex')
ylabel('$f_{\rm{\theta}}$','Interpreter','latex')



subplot(1,2,2);
X = minKeq_vals;
histogram(X,1000,'Normalization','pdf','FaceColor','c','FaceAlpha',1,'EdgeColor','c');hold on;
xlabel('$\rm{K}$','Interpreter','latex')
ylabel('$f_{\rm{K}}$','Interpreter','latex')




%%

clear all; clc
N = 1e7;  K0 = 1;
planes = [1  1  0;
          1 -1  0;
         -1  1  0;
         -1 -1  0;
          1  0  1;
          1  0 -1;
         -1  0  1;
         -1  0 -1;
          0  1  1;
          0  1 -1;
          0 -1  1;
          0 -1 -1]'/sqrt(2);
minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C = abs(R(3,:)*planes);
    K_eq = K0./max(eps,C);
    minKeq_vals(i) = min(K_eq);
end

% figure;
% X = acos(1./minKeq_vals);
% histogram(X,1000,'Normalization','pdf','FaceColor','y','FaceAlpha',1,'EdgeColor','y');hold on;
% 
% 
% alpha=linspace(0,pi/4,1000);
% phi0=asin(1/sqrt(3));
% N1=1/(1 - pi/(4*acos(1/sqrt(3))));
% N2=1/(1 - pi/(6*asin(1/sqrt(3))));
% 
% f=zeros(size(alpha));
% for i=1:length(alpha)
%     a=alpha(i);
%     if a<=pi/6
%         f(i)=((pi-2*phi0)*N1/pi + 2*phi0*N2/pi)*sin(a);
%     elseif a<=phi0
%         f(i)=(N1/pi*(pi - 2*acos((cos(a)/sin(a))/sqrt(3)) - 2*phi0) + ...
%               N2/pi*(2*phi0 - 2*acos((cos(a)/sin(a))/sqrt(3))))*sin(a);
%     else
%         f(i)= (N1/pi*(pi - 2*acos((cos(a)/sin(a))/sqrt(3)) - 2*phi0) )*sin(a);
%     end
% end
% 
% plot(alpha,f,'r','LineWidth',2); hold on;
% 
% xlabel('\theta')
% ylabel('PDF')


figure;
subplot(1,2,1);
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','y','FaceAlpha',1,'EdgeColor','b');hold on;


alpha=linspace(0,pi/4,1000);
phi0=asin(1/sqrt(3));
N1=1/(1 - pi/(4*acos(1/sqrt(3))));
N2=1/(1 - pi/(6*asin(1/sqrt(3))));

f=zeros(size(alpha));
for i=1:length(alpha)
    a=alpha(i);
    if a<=pi/6
        f(i)=((pi-2*phi0)*N1/pi + 2*phi0*N2/pi)*sin(a);
    elseif a<=phi0
        f(i)=(N1/pi*(pi - 2*acos((cos(a)/sin(a))/sqrt(3)) - 2*phi0) + ...
              N2/pi*(2*phi0 - 2*acos((cos(a)/sin(a))/sqrt(3))))*sin(a);
    else
        f(i)= (N1/pi*(pi - 2*acos((cos(a)/sin(a))/sqrt(3)) - 2*phi0) )*sin(a);
    end
end

plot(alpha,f,'r','LineWidth',2); hold on;

xlabel('$\theta$','Interpreter','latex')
ylabel('$f_{\rm{\theta}}$','Interpreter','latex')



subplot(1,2,2);
X = minKeq_vals;
histogram(X,1000,'Normalization','pdf','FaceColor','y','FaceAlpha',1,'EdgeColor','b');hold on;
xlabel('$\rm{K}$','Interpreter','latex')
ylabel('$f_{\rm{K}}$','Interpreter','latex')


%%

clear all; clc
N = 1e7;  K0 = 1;
planes = [ 1  1  1;
          -1  1  1;
           1 -1  1;
           1  1 -1;
          -1 -1  1;
          -1  1 -1;
           1 -1 -1;
          -1 -1 -1 ]'/sqrt(3);

minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C = abs(R(3,:)*planes);
    K_eq = K0./max(eps,C);
    minKeq_vals(i) = min(K_eq);
end

figure;
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','c','FaceAlpha',1,'EdgeColor','c');hold on;



%%

clear all; clc
N = 1e7;  K0 = 1;
planes = [1  1  0;
          1 -1  0;
         -1  1  0;
         -1 -1  0;
          1  0  1;
          1  0 -1;
         -1  0  1;
         -1  0 -1;
          0  1  1;
          0  1 -1;
          0 -1  1;
          0 -1 -1]'/sqrt(2);

planes1 = [1  0  0;
          -1  0  0;
           0  1  0;
           0 -1  0;
           0  0  1;
           0  0 -1]';





planes = [planes  planes1];

minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C = abs(R(3,:)*planes);
    K_eq = K0./max(eps,C);
    minKeq_vals(i) = min(K_eq);
end

figure;
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','c','FaceAlpha',1,'EdgeColor','c');hold on;


%%



clear all; clc
N = 1e7;  K0 = 1;
planes = [ 1  1  1;
          -1  1  1;
           1 -1  1;
           1  1 -1;
          -1 -1  1;
          -1  1 -1;
           1 -1 -1;
          -1 -1 -1 ]'/sqrt(3);

planes1 = [1  0  0;
          -1  0  0;
           0  1  0;
           0 -1  0;
           0  0  1;
           0  0 -1]';
planes = [planes  planes1];

minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C = abs(R(3,:)*planes);
    K_eq = K0./max(eps,C);
    minKeq_vals(i) = min(K_eq);
end

figure;
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','c','FaceAlpha',1,'EdgeColor','c');hold on;

%%


clear all; clc
N = 1e7;  K0 = 1;
planes = [ 1  1  1;
          -1  1  1;
           1 -1  1;
           1  1 -1;
          -1 -1  1;
          -1  1 -1;
           1 -1 -1;
          -1 -1 -1 ]'/sqrt(3);

planes1 = [1  0  0;
          -1  0  0;
           0  1  0;
           0 -1  0;
           0  0  1;
           0  0 -1]';


planes2 = [1  1  0;
          1 -1  0;
         -1  1  0;
         -1 -1  0;
          1  0  1;
          1  0 -1;
         -1  0  1;
         -1  0 -1;
          0  1  1;
          0  1 -1;
          0 -1  1;
          0 -1 -1]'/sqrt(2);


planes = [planes  planes1  planes2];

minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C = abs(R(3,:)*planes);
    K_eq = K0./max(eps,C);
    minKeq_vals(i) = min(K_eq);
end

figure;
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','c','FaceAlpha',1,'EdgeColor','c');hold on;



%%

clear all; clc;
N = 1e6;  R1 = sqrt(2); R2 = sqrt(1);

planes1 = [1  0  0;
          -1  0  0;
           0  1  0;
           0 -1  0;
           0  0  1;
           0  0 -1]';

planes2 = [1  1  0;
          1 -1  0;
         -1  1  0;
         -1 -1  0;
          1  0  1;
          1  0 -1;
         -1  0  1;
         -1  0 -1;
          0  1  1;
          0  1 -1;
          0 -1  1;
          0 -1 -1]'/sqrt(2);

planes = [planes1  planes2];
minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C1 = abs(R(3,:)*planes1);
    C2 = abs(R(3,:)*planes2);

    K1 = R1./max(eps,C1);
    K2 = R2./max(eps,C2);
    K_eq = [K1 K2];

    minKeq_vals(i) = min(K_eq);
end

figure;
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','r','FaceAlpha',1,'EdgeColor','r');hold on;



%%



clear all; clc;


N = 1e6;  R1 = sqrt(X); R2 = sqrt(1);

planes1 = [1  0  0;
          -1  0  0;
           0  1  0;
           0 -1  0;
           0  0  1;
           0  0 -1]';


planes2 = [1  1  0;
          1 -1  0;
         -1  1  0;
         -1 -1  0;
          1  0  1;
          1  0 -1;
         -1  0  1;
         -1  0 -1;
          0  1  1;
          0  1 -1;
          0 -1  1;
          0 -1 -1]'/sqrt(2);


planes = [planes1  planes2];

minKeq_vals = zeros(N,1);
z_hat = [0;0;1];

for i = 1:N
    q = randn(4,1); q = q/norm(q);
    w = q(1); x_ = q(2); y_ = q(3); z_ = q(4);
    R = [1-2*(y_^2+z_^2),   2*(x_*y_-w*z_),    2*(x_*z_+w*y_);
         2*(x_*y_+w*z_),   1-2*(x_^2+z_^2),    2*(y_*z_-w*x_);
         2*(x_*z_-w*y_),    2*(y_*z_+w*x_),   1-2*(x_^2+y_^2)];
    C1 = abs(R(3,:)*planes1);
    C2 = abs(R(3,:)*planes2);

    K1 = R1./max(eps,C1);
    K2 = R2./max(eps,C2);
    K_eq = [K1 K2];

    minKeq_vals(i) = min(K_eq);
end

figure;
X = acos(1./minKeq_vals);
histogram(X,1000,'Normalization','pdf','FaceColor','r','FaceAlpha',1,'EdgeColor','r');hold on;




%%


clear; clc;
planes1 = [ 1  0  0;
           -1  0  0;
            0  1  0;
            0 -1  0;
            0  0  1;
            0  0 -1 ]';

planes2 = [ 1  1  0;
            1 -1  0;
           -1  1  0;
           -1 -1  0;
            1  0  1;
            1  0 -1;
           -1  0  1;
           -1  0 -1;
            0  1  1;
            0  1 -1;
            0 -1  1;
            0 -1 -1 ]' / sqrt(2);

N  = 1e7;                      
Ps = 1 : 0.1 : 2.15;               
nP = numel(Ps);

nRows = 3;
nCols = ceil(nP / nRows);
figure('Color','w');

for k = 1:nP
    Pval = Ps(k);               
    R1   = sqrt(Pval);          % ⟨100⟩ coefficient
    R2   = 1;                   % ⟨110⟩ coefficient

    minKeq_vals = zeros(N,1);

    for i = 1:N
        q = randn(4,1);  q = q / norm(q);     
        w = q(1);  x = q(2);  y = q(3);  z = q(4);

        R = [ 1-2*(y^2+z^2),    2*(x*y-w*z),    2*(x*z+w*y);
              2*(x*y+w*z),    1-2*(x^2+z^2),    2*(y*z-w*x);
              2*(x*z-w*y),      2*(y*z+w*x),  1-2*(x^2+y^2) ];

        C1 = abs(R(3,:) * planes1);      % 1×6
        C2 = abs(R(3,:) * planes2);      % 1×12

        K1 = R1 ./ C1;                   % 1×6
        K2 = R2 ./ C2;                   % 1×12
        minKeq_vals(i) = min( [K1 K2] );
    end

    subplot(nRows,nCols,k)
    theta = acos(1 ./ minKeq_vals);      
    histogram(theta, 1000, ...
              'Normalization','pdf', ...
              'EdgeColor','none', ...
              'FaceColor',[0.85 0.1 0.1], ...
              'FaceAlpha',0.6);
    title(sprintf('P = %.1f', Pval))
    xlabel('$\theta_{\min}$','Interpreter','latex')
    ylabel('pdf')
end


%%


clear; clc;
planes1 = [ 1  0  0;
           -1  0  0;
            0  1  0;
            0 -1  0;
            0  0  1;
            0  0 -1 ]';

planes2 = [ 1  1  0;
            1 -1  0;
           -1  1  0;
           -1 -1  0;
            1  0  1;
            1  0 -1;
           -1  0  1;
           -1  0 -1;
            0  1  1;
            0  1 -1;
            0 -1  1;
            0 -1 -1 ]' / sqrt(2);

N  = 1e7;                      
Ps = 1 : 0.1 : 2.15;               
nP = numel(Ps);

nRows = 3;
nCols = ceil(nP / nRows);
figure('Color','w');

for k = 1:nP
    Pval = Ps(k);               
    R2   = sqrt(Pval);          % ⟨100⟩ coefficient
    R1   = 1;                   % ⟨110⟩ coefficient

    minKeq_vals = zeros(N,1);

    for i = 1:N
        q = randn(4,1);  q = q / norm(q);     
        w = q(1);  x = q(2);  y = q(3);  z = q(4);

        R = [ 1-2*(y^2+z^2),    2*(x*y-w*z),    2*(x*z+w*y);
              2*(x*y+w*z),    1-2*(x^2+z^2),    2*(y*z-w*x);
              2*(x*z-w*y),      2*(y*z+w*x),  1-2*(x^2+y^2) ];

        C1 = abs(R(3,:) * planes1);      % 1×6
        C2 = abs(R(3,:) * planes2);      % 1×12

        K1 = R1 ./ C1;                   % 1×6
        K2 = R2 ./ C2;                   % 1×12
        minKeq_vals(i) = min( [K1 K2] );
    end

    subplot(nRows,nCols,k)
    theta = acos(1 ./ minKeq_vals);      
    histogram(theta, 1000, ...
              'Normalization','pdf', ...
              'EdgeColor','none', ...
              'FaceColor','c', ...
              'FaceAlpha',0.6);
    title(sprintf('Q = %.1f', Pval))
    xlabel('$\theta_{\min}$','Interpreter','latex')
    ylabel('pdf')
   % xlim([0 pi/4])
end


%%


clear; clc
planes1=[ 1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]';
planes2=[ 1 1 0;1 -1 0;-1 1 0;-1 -1 0;1 0 1;1 0 -1;-1 0 1;-1 0 -1;0 1 1;0 1 -1;0 -1 1;0 -1 -1]'/sqrt(2);
N=1e5;
ratios=linspace(2,0.5,9);
figure('Color','w')
for k=1:length(ratios)
    r=ratios(k);
    if r>=1
        R100=sqrt(r); R110=1;
    else
        R100=1; R110=sqrt(1/r);
    end
    Kmin=zeros(N,1);
    for i=1:N
        q=randn(4,1); q=q/norm(q);
        w=q(1); x=q(2); y=q(3); z=q(4);
        R=[1-2*(y^2+z^2) 2*(x*y-w*z) 2*(x*z+w*y);
           2*(x*y+w*z) 1-2*(x^2+z^2) 2*(y*z-w*x);
           2*(x*z-w*y) 2*(y*z+w*x) 1-2*(x^2+y^2)];
        C1=abs(R(3,:)*planes1);
        C2=abs(R(3,:)*planes2);
        Kmin(i)=min([R100./C1 R110./C2]);
        
    end
    subplot(3,3,k)
    histogram(acos(1./Kmin),1000,'Normalization','pdf','EdgeColor','none')
    title(sprintf('P/Q = %.2f',r))
    xlabel('$\theta_{\min}$','Interpreter','latex')
    ylabel('$f$','Interpreter','latex')
end


%%


clear; clc
planes1=[1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]';
planes2=[1 1 0; 1 -1 0; -1 1 0; -1 -1 0; 1 0 1; 1 0 -1; -1 0 1; -1 0 -1; 0 1 1; 0 1 -1; 0 -1 1; 0 -1 -1]'/sqrt(2);
N=1e5;
ratios=linspace(2,0.5,9);
figure('Color','w')
for k=1:length(ratios)
    r=ratios(k);
    if r>=1, R100=sqrt(r); R110=1; else, R100=1; R110=sqrt(1/r); end

    Kmin_vals = zeros(N,1);
    K11 = zeros(N,1);
    K22 = zeros(N,1);
    for i = 1:N
        q = randn(4,1); q = q/norm(q);
        w = q(1); x = q(2); y = q(3); z = q(4);
        R = [1-2*(y^2+z^2)  2*(x*y-w*z)  2*(x*z+w*y);
            2*(x*y+w*z)  1-2*(x^2+z^2)  2*(y*z-w*x);
            2*(x*z-w*y)    2*(y*z+w*x)  1-2*(x^2+y^2)];
        C1 = abs(R(3,:)*planes1);
        C2 = abs(R(3,:)*planes2);
        K1 = R100./C1;
        K2 = R110./C2;
        [Kmin, idx] = min([K1, K2]);
        Kmin_vals(i) = Kmin;
        if idx <= numel(K1)
            K11(i) = Kmin;
        else
            K22(i) = Kmin;
        end
    end


    subplot(3,3,k)
    histogram(acos(1./Kmin_vals),1000,'Normalization','pdf','FaceColor','b','EdgeColor','none')
    title(sprintf('P/Q = %.2f',r))
    xlabel('$\theta_{\min}$','Interpreter','latex')
    ylabel('$f$','Interpreter','latex')
end


%%

clear; clc
planes1=[1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]';
planes2=[1 1 0; 1 -1 0; -1 1 0; -1 -1 0; 1 0 1; 1 0 -1; -1 0 1; -1 0 -1; 0 1 1; 0 1 -1; 0 -1 1; 0 -1 -1]'/sqrt(2);
N=1e7;
ratios=linspace(2,0.5,9);

figure('Color','w');
hold on;
for k=1:length(ratios)
    r=ratios(k);
    if r>=1, R100=sqrt(r); R110=1; else, R100=1; R110=sqrt(1/r); end
    Kmin=zeros(N,1);
    K11=zeros(N,1);
    K22=zeros(N,1);
    for i=1:N
        q=randn(4,1); q=q/norm(q);
        w=q(1); x=q(2); y=q(3); z=q(4);
        R=[1-2*(y^2+z^2) 2*(x*y-w*z) 2*(x*z+w*y);
           2*(x*y+w*z) 1-2*(x^2+z^2) 2*(y*z-w*x);
           2*(x*z-w*y) 2*(y*z+w*x) 1-2*(x^2+y^2)];
        C1=abs(R(3,:)*planes1);
        C2=abs(R(3,:)*planes2);
        K1=R100./C1;
        K2=R110./C2;
        [m,idx]=min([K1,K2]);
        Kmin(i)=m;
        [KmiN, ~] = min([K1, K2]);
        Kmin_vals(i) = KmiN;
        if idx<=6
            K11(i)=m;
        else
            K22(i)=m;
        end
    end
    subplot(3,3,k); hold on;
     % histogram(acos(1./K11(K11>0)),1000,'Normalization','countdensity','FaceColor','c','FaceAlpha',0.9,'EdgeColor','none'); 
     % histogram(acos(1./K22(K22>0)),1000,'Normalization','countdensity','FaceColor','r','FaceAlpha',0.35,'EdgeColor','none')
     % histogram(acos(1./Kmin_vals),1000,'Normalization','countdensity','FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.7,'EdgeColor','none'); hold on;
     % xlim([0 1])

     % histogram(acos(1./K11(K11>0)),1000,'Normalization','countdensity','FaceColor','c','FaceAlpha',0.9,'EdgeColor','none'); 
     % histogram(acos(1./K22(K22>0)),1000,'Normalization','countdensity','FaceColor','r','FaceAlpha',0.35,'EdgeColor','none')
    histogram(Kmin_vals,1000,'Normalization','countdensity','FaceColor','k','FaceAlpha',0.7,'EdgeColor','none'); hold on;
    histogram(K11(K11>0),1000,'Normalization','countdensity','FaceColor','c','FaceAlpha',0.5,'EdgeColor','none'); hold on;
    histogram(K22(K22>0),1000,'Normalization','countdensity','FaceColor','r','FaceAlpha',0.3,'EdgeColor','none'); hold on;

    hold on
    title(sprintf('P/Q = %.2f',r))
    xlabel('K')
    ylabel('Frequency')
    xlim([1 max(Kmin_vals)])
    drawnow;

end

%%

clear; clc
planes1=[1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]';
planes2=[1 1 0; 1 -1 0; -1 1 0; -1 -1 0; 1 0 1; 1 0 -1; -1 0 1; -1 0 -1; 0 1 1; 0 1 -1; 0 -1 1; 0 -1 -1]'/sqrt(2);
N=1e7;
ratios=linspace(2,0.5,20);

figure('Color','w');
hold on;
for k=1:length(ratios)
    r=ratios(k);
    if r>=1, R100=sqrt(r); R110=1; else, R100=1; R110=sqrt(1/r); end
    Kmin=zeros(N,1);
    K11=zeros(N,1);
    K22=zeros(N,1);
    for i=1:N
        q=randn(4,1); q=q/norm(q);
        w=q(1); x=q(2); y=q(3); z=q(4);
        R=[1-2*(y^2+z^2) 2*(x*y-w*z) 2*(x*z+w*y);
           2*(x*y+w*z) 1-2*(x^2+z^2) 2*(y*z-w*x);
           2*(x*z-w*y) 2*(y*z+w*x) 1-2*(x^2+y^2)];
        C1=abs(R(3,:)*planes1);
        C2=abs(R(3,:)*planes2);
        K1=R100./C1;
        K2=R110./C2;
        [m,idx]=min([K1,K2]);
        Kmin(i)=m;
        [KmiN, ~] = min([K1, K2]);
        Kmin_vals(i) = KmiN;
        if idx<=6
            K11(i)=m;
        else
            K22(i)=m;
        end
    end
     disp(length(K11(K11>0))/N);
     disp(length(K22(K22>0))/N);
     disp('______');

     plot(r, length(K11(K11>0))/N,'rs'); hold on;
     plot(r, length(K22(K22>0))/N,'bs'); hold on;

    % hold on
    % title(sprintf('P/Q = %.2f',r))
     xlabel('P/Q')
     ylabel('Fraction')
    % xlim([1 max(Kmin_vals)])
    % drawnow;

end
