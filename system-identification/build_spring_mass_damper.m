function [A,B,C,D] = build_spring_mass_damper(par, T, aux)
% par: [k, b, qnoise]
% aux: [dt, m, rnoise]
stiff = par(1);
damp = par(2);
mass = par(3);
% X0 = [par(3); par(4)];
% qnoise = abs(par(3));

dt = aux(1);
% m = aux(2);
% R2 = dt * abs(aux(3)); % Known measurement noise variance

% spring-mass-damper
% A = expm([0, 1; -stiff/mass, -damp/mass]*dt);
% B = [0; 1/mass];  % no input
% C = [1, 0];  % measure displacement
% D = 0;

A = [0, 1; -stiff/mass, -damp/mass];
B = [0; 1/mass];  % no input
C = [1, 0];  % measure displacement
D = 0;

% R1 = [  dt^3/3, dt^2/2
%         dt^2/2, dt] * qnoise;  % acceleration noise
% nx = 2;
% ny = 1;

% sys = ss(A, eye(nx), C, zeros(ny,0), dt);
% R12 = zeros(nx, ny);
% [~,K] = kalman(sys, R1, R2, R12);





% Fs = 1/dt;
% forceCov = 0.5e-1 ^ 2;  % force noise covariance
% mocapCov = 2e-3^2;  % mocap noise covariance
% uCov = diag([mocapCov, forceCov]);
% T = diag([1e3, 1e3, 1e3]); % eye(4); %diag([1, 1/Fs, 1/Fs^2]);  % state transformation
% A = T * [0, 0, 0
%          1, 0, 0
%          0, 1, 0] / T;
% B = T * [1, 0
%          0, 0
%          0, 0];
% C = -[stiff, damp*Fs, 1*mass*Fs^2] * [[   0,  1,    0]
%                                     [ 1/2,  0, -1/2]
%                                     [   1, -2,    1]] / T;
% D = [0, 1];
% Q = B * uCov * B'; %zeros(3);
% R = 1e-3 ^ 2 + D * uCov * D';  % trust on the equation of motion, given in Nm



% Fs = 1/dt;
% forceCov = 0.5e-1 ^ 2;  % force noise covariance
% mocapCov = 2e-3^2;  % mocap noise covariance
% uCov = diag([mocapCov]);
% T = diag([1e3, 1e3, 1e3]); % eye(4); %diag([1, 1/Fs, 1/Fs^2]);  % state transformation
% A = T * [0, 0, 0
%          1, 0, 0
%          0, 1, 0] / T;
% B = T * [1
%          0
%          0];
% C = -[stiff, damp*Fs, 1*mass*Fs^2] * [[   0,  1,    0]
%                                     [ 1/2,  0, -1/2]
%                                     [   1, -2,    1]] / T;
% D = [0];
% Q = B * uCov * B'; %zeros(3);
% R = 1e-3 ^ 2 + D * uCov * D' + forceCov;  % trust on the equation of motion, given in Nm
