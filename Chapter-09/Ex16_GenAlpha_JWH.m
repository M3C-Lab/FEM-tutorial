% clean the memory and the screen
clear all; clc; close all;

% Setup the alpha parameters
rho_inf = 0.5;
alpha_m = 0.5 * (3.0 - rho_inf) / (1.0 + rho_inf);
alpha_f = 1.0 / (1.0 + rho_inf);
gamma   = 0.5 - alpha_f + alpha_m;

% Setup the physical problem parameters
% mass
m1 = 1.0;
m2 = 1.0;

% stiffness
k1 = 1.0e4;
k2 = 1.0;

% setup the matrices
M = [ m1, 0.0; 0.0, m2 ];

K = [ (k1+k2), -k2; -k2, k2 ];

% setup the initial condition
d0 = [ 1 ; 10 ];

v0 = [ 0; 0 ];

% determine the natural frequencies
lambda = eig(K, M);
omega = sqrt(lambda);

T1 = 2 * pi / omega(1);
T2 = 2 * pi / omega(2);

dt = T1 / 200;

T_final = 5 * T1;

N = ceil(T_final / dt);

% initial acceleration
a0 = M \ (-K * d0);

% allocate solutions
v = zeros(2, N+1); d = v; dot_d = v; dot_v = v;

v(:,1) = v0; dot_v(:,1) = a0;
d(:,1) = d0; dot_d(:,1) = v0;

% Matrix
LEFT = alpha_m * M + ( alpha_f^2 * gamma^2 * dt^2 / alpha_m ) * K;

for n = 2 : N+1
  % predictor
  d_tilde = d(:,n-1) + (1-gamma) * dt * dot_d(:, n-1);
  v_tilde = v(:,n-1) + (1-gamma) * dt * dot_v(:, n-1);
  
  RIGHT = (1-alpha_m) * M * dot_v(:,n-1) + (1-alpha_f) * K * d(:,n-1) ...
    + alpha_f * K * d_tilde + (alpha_f * gamma * dt / alpha_m) ...
    * K * ( (1-alpha_f)*v(:,n-1) - (1-alpha_m) * dot_d(:,n-1) + alpha_f * v_tilde );

  RIGHT = -RIGHT;

  dot_v(:,n) = LEFT \ RIGHT;

  v(:,n) = v_tilde + gamma * dt * dot_v(:,n);

  dot_d(:,n) = ( (1-alpha_f)*v(:,n-1) + alpha_f*v(:,n) - (1-alpha_m)*dot_d(:, n-1) ) / alpha_m;

  d(:,n) = d_tilde + gamma * dt * dot_d(:,n);
end

% visualization
t = 1: 1 : N;

subplot(2,2,1), plot(t, d(1,1:N)); grid on;
subplot(2,2,2), plot(t, d(2,1:N)); grid on;
subplot(2,2,3), plot(t, v(1,1:N)); grid on;
subplot(2,2,4), plot(t, v(2,1:N)); grid on;

% EOF