function [A,B,C,D,K,X0] = mynoise(par,T,aux)

dt = aux(2);
% R2 = dt*diag([1e-4, 0.9e-4, 1.1e-4]); % Known measurement noise variance
R2 = dt*diag(1e-4); % Known measurement noise variance

% spring-mass-damper
A = expm([0, 1; -100/0.1, -3/0.1]*dt);
B = zeros(2, 0);
C = [1, 0];
D = zeros(1, 0);
R1 = [   dt^3/3, dt^2/2
        dt^2/2, dt] * aux(1);
nx = 2;
ny = 1;

% % constant accel
% A = [1, dt, dt^2/2
%     0,  1,  dt
%     0,  0,  1];
% B = zeros(3, 0);
% C = [1, 0, 0; 0.5, 0, 0; -1.2, 0, 0];
% D = zeros(3, 0);
% R1 = [   dt^5/20,   dt^4/8,  dt^3/6
%         dt^4/8,     dt^3/3, dt^2/2
%         dt^3/6,     dt^2/2, dt] * par(1);  % white noise in accel

sys = ss(A, eye(nx), C, zeros(ny,0), dt);
R12 = zeros(1, 1);
[~,K] = kalman(sys, R1, R2, R12);
X0 = [par(2); par(3)];
