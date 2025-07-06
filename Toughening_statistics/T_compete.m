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


