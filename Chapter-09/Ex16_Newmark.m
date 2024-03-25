% clean the memory and the screen
clear all; clc; close all;

% Setup the Newmark parameters
% Central difference: beta = 0, gamma = 0.5
% Trapezoidal: beta = 0.25, gamma = 0.5
% Damped Newmark: beta = 0.3025, gamma = 0.6
beta = 0.0; 
gamma = 0.5;

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

dt = T1 / 20;

T_final = 5 * T1;

N = ceil(T_final / dt);

% initial acceleration
a0 = M \ (-K * d0);

% allocate solutions
a = zeros(2, N+1);
v = a;
d = a;

a(:,1) = a0; v(:,1) = v0; d(:,1) = d0;

% matrix
LEFT = M + beta * dt * dt * K;

for n = 2 : N+1
    % predictor
    d_tilde = d(:,n-1) + dt * v(:,n-1) + 0.5 * dt * dt * (1 - 2 * beta) * a(:, n-1);
    v_tilde = v(:,n-1) + (1-gamma) * dt * a(:, n-1);
    
    a(:,n) = LEFT \ (-K * d_tilde);
    
    d(:,n) = d_tilde + beta * dt * dt * a(:,n);
    v(:,n) = v_tilde + gamma * dt * a(:,n);
end

% visualization
t = 1: 1 : N;

subplot(2,2,1), plot(t, d(1,1:N)); grid on;
subplot(2,2,2), plot(t, d(2,1:N)); grid on;
subplot(2,2,3), plot(t, v(1,1:N)); grid on;
subplot(2,2,4), plot(t, v(2,1:N)); grid on;

% EOF