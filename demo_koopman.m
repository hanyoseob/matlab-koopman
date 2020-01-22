clear all;
close all;
clc;

%% generate grid geometry
% space
L = 30;
n = 512;
x = linspace(-L/2,L/2,n+1);
x = x(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1].';

% time
slices = 20;
t = linspace(0, pi, slices + 1);
dt = t(2) - t(1);

nf = 0.5;
slicesf = 2*slices*nf;
tf = linspace(0, 2*pi*nf, slicesf + 1);

% initial conditions
N = 2;
u0 = N*(sech(x)).';
ut0 = fft(u0);

%% create two spatio-temoporal data by solving ODE solver (ode45)
[t, utsol] = ode45('dmd_soliton_rhs', t, ut0, [], k);
usol = ifft(utsol, [], 2);

%% dynamic mode decomposition (DMD)
X = usol.';  % in C^(spatio, temporal)
r = 10;

[X_dmd, Phi, omega] = dmd_func(X, t, dt, r, u0);

% ERROR 
for j=1:length(t); error_dmd(j)=norm(X_dmd(:,j)-usol(j,:).'); end

%% Koopman operator; correct data extention
Y1 = [X; (X.*abs(X).^2)];
u0_1 = [u0; (u0.*abs(u0).^2)];
r = 10;

[Y_koop1, Phi_y1, omega_x1] = dmd_func(Y1, t, dt, r, u0_1);
X_koop1 = Y_koop1(1:n, :);
Phi_x1 = Phi_y1(1:n, :);

% ERROR 
for j=1:length(t); error_koop1(j)=norm(X_koop1(:,j)-usol(j,:).'); end

%% Koopman operator; uncorrect data extention
Y2 = [X; (abs(X).^2)];
u0_2 = [u0; (abs(u0).^2)];
r = 10;

[Y_koop2, Phi_y2, omega_x2] = dmd_func(Y2, t, dt, r, u0_2);
X_koop2 = Y_koop2(1:n, :);
Phi_x2 = Phi_y2(1:n, :);

% ERROR 
for j=1:length(t); error_koop2(j)=norm(X_koop2(:,j)-usol(j,:).'); end

%% Display reconstruction results
figure('name', 'Reconstruction results');
subplot(2,2,1); waterfall(x,t,abs(usol)); colormap([0, 0 ,0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('time', 'Fontsize', [10]); title('Measurements');
set(gca, 'Fontsize', [10], 'Ylim', [0, pi], 'Ytick', [0, 3, 6], 'Xlim', [-15, 15], 'Xtick',[ -15, 0, 15], 'Zlim', [0, 4], 'Ztick', [0, 2, 4]);

subplot(2,2,2); waterfall(x,tf,abs(X_dmd).'); colormap([0, 0, 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('time', 'Fontsize', [10]); title('DMD reconstructions');
set(gca, 'Fontsize', [10], 'Ylim', [0, pi], 'Ytick', [0, 3, 6], 'Xlim', [-15, 15], 'Xtick',[ -15, 0, 15], 'Zlim', [0, 4], 'Ztick', [0, 2, 4]);

subplot(2,2,3); waterfall(x,tf,abs(X_koop1).'); colormap([0, 0, 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('time', 'Fontsize', [10]); title({'Koopman reconstructions', 'Good Embedding'});
set(gca, 'Fontsize', [10], 'Ylim', [0, pi], 'Ytick', [0, 3, 6], 'Xlim', [-15, 15], 'Xtick',[ -15, 0, 15], 'Zlim', [0, 4], 'Ztick', [0, 2, 4]);

subplot(2,2,4); waterfall(x,tf,abs(X_koop2).'); colormap([0, 0, 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('time', 'Fontsize', [10]); title({'Koopman reconstructions', 'Bad Embedding'});
set(gca, 'Fontsize', [10], 'Ylim', [0, pi], 'Ytick', [0, 3, 6], 'Xlim', [-15, 15], 'Xtick',[ -15, 0, 15], 'Zlim', [0, 4], 'Ztick', [0, 2, 4]);

%% Display Eigen values
figure('name', 'Eigen values');
subplot(3,3,1), waterfall(x,tf,abs(X_dmd).'); colormap([0, 0, 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('time', 'Fontsize', [10]); title('DMD reconstructions');
set(gca, 'Fontsize', [10], 'Ylim', [0, pi], 'Ytick', [0, 3, 6], 'Xlim', [-15, 15], 'Xtick',[ -15, 0, 15], 'Zlim', [0, 4], 'Ztick', [0, 2, 4]);

subplot(3,3,2), waterfall(x,tf,abs(X_koop1).'); colormap([0, 0, 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('time', 'Fontsize', [10]); title({'Koopman reconstructions', 'Good Embedding'});
set(gca, 'Fontsize', [10], 'Ylim', [0, pi], 'Ytick', [0, 3, 6], 'Xlim', [-15, 15], 'Xtick',[ -15, 0, 15], 'Zlim', [0, 4], 'Ztick', [0, 2, 4]);

subplot(3,3,3), waterfall(x,tf,abs(X_koop2).'); colormap([0, 0, 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('time', 'Fontsize', [10]); title({'Koopman reconstructions', 'Bad Embedding'});
set(gca, 'Fontsize', [10], 'Ylim', [0, pi], 'Ytick', [0, 3, 6], 'Xlim', [-15, 15], 'Xtick',[ -15, 0, 15], 'Zlim', [0, 4], 'Ztick', [0, 2, 4]);

subplot(3,3,4); plot(omega, 'ko', 'Linewidth', [2]); grid on; axis([-5, 1, -20, 20]);
xlabel('Re', 'Fontsize', [10]); ylabel('Im', 'Fontsize', [10]); title('Eigen value_{DMD}');
set(gca, 'Fontsize', [10]);

subplot(3,3,5); plot(omega_x1, 'ko', 'Linewidth', [2]); grid on; axis([-5, 1, -20, 20]);
xlabel('Re', 'Fontsize', [10]); ylabel('Im', 'Fontsize', [10]); title('Eigen value_{Koopman; Good}');
set(gca, 'Fontsize', [10]);

subplot(3,3,6); plot(omega_x2, 'ko', 'Linewidth', [2]); grid on; axis([-5, 1, -20, 20]);
xlabel('Re', 'Fontsize', [10]); ylabel('Im', 'Fontsize', [10]); title('Eigen value_{Koopman; Bad}');
set(gca, 'Fontsize', [10]);

subplot(3,3,[7,9]); semilogy(t,error_dmd,'k-',t,error_koop1,'k:',t,error_koop2,'k--','Linewidth',[2]);
xlabel('time', 'Fontsize', [10]); ylabel('Error', 'Fontsize', [10]); title('Mean Square Error');
set(gca,'Fontsize', [10], 'Ylim', [10^(-6), 10^2], 'Ytick', [10^(-6), 10^(-4), 10^(-2), 10^0, 10^2]);
axis([0, pi, 10^(-6), 10^2]); grid on; grid minor;
legend('DMD', 'koopman; good', 'koopman; bad');

%% Display Eigen vectors
figure('name', 'Eigen vectors');
subplot(1,3,1); waterfall(x,1:r,abs(Phi).'); colormap([0 0 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('rank', 'Fontsize', [10]); title('Eigen vector_{DMD}')
set(gca, 'Fontsize', [10], 'Ylim', [1, r], 'Ytick', [1, r], 'Xlim', [-8, 8], 'Xtick',[ -8, 0, 8], 'Zlim', [0 0.4], 'Ztick', [0 0.2 0.4]);

subplot(1,3,2); waterfall(x,1:r,abs(Phi_x1).'); colormap([0 0 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('rank', 'Fontsize', [10]); title('Eigen vector_{Koopman; Good}')
set(gca, 'Fontsize', [10], 'Ylim', [1, r], 'Ytick', [1, r], 'Xlim', [-8, 8], 'Xtick',[ -8, 0, 8], 'Zlim', [0 0.05], 'Ztick', [0 0.025 0.05]);

subplot(1,3,3); waterfall(x,1:r,abs(Phi_x2).'); colormap([0 0 0]);
xlabel('\xi', 'Fontsize', [10]); ylabel('rank', 'Fontsize', [10]); title('Eigen vector_{Koopman; Bad}')
set(gca, 'Fontsize', [10], 'Ylim', [1, r], 'Ytick', [1, r], 'Xlim', [-8, 8], 'Xtick',[ -8, 0, 8], 'Zlim', [0 0.3], 'Ztick', [0 0.15 0.3]);

