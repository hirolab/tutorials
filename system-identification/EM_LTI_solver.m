function cost = EM_LTI_solver(params, xv, fv)

% Tp = eye(3); %diag([1/1e2, 1/1e0, 1/1e-1]);
% params = Tp \ params;

stiff = params(1);
damp = params(2);
mass = params(3);

forceCov = 0.5e-1 ^ 2;  % force noise covariance
mocapCov = 2e-3^2;  % mocap noise covariance
Fs = 183;

% impedance D =======================
u = [xv(3:end), fv(1:end-2)];
uCov = diag([mocapCov, forceCov]);
y = zeros(size(u,1),1);
T = diag([1e3, 1e3, 1e3]); % eye(4); %diag([1, 1/Fs, 1/Fs^2]);  % state transformation
A = T * [0, 0, 0
         1, 0, 0
         0, 1, 0] / T;
B = T * [1, 0
         0, 0
         0, 0];
C = -[stiff, damp*Fs, 1*mass*Fs^2] * [[   0,  1,    0]
                                    [ 1/2,  0, -1/2]
                                    [   1, -2,    1]] / T;
D = [0, 1];
Q = B * uCov * B'; %zeros(3);
R = 1e-3 ^ 2 + D * uCov * D';  % trust on the equation of motion, given in Nm
xk = T * [xv(1); xv(1); xv(1)];
P = T * diag([mocapCov, mocapCov, mocapCov]) * T';

% admittance =======================
% u = 0.5*(fv + fv([2:end, end]));  % shift force 0.5 sample
% uCov = forceCov * 2 * 0.5^2;
% y = xv;
% T = diag([1e3, 1e3/Fs]);
% sys = c2d(ss([0, 1; -stiff/mass, -damp/mass], [0; 1/mass], [1, 0], 0), 1/Fs);
% A = T*sys.A/T; B = T*sys.B; C = sys.C/T; D = sys.D;
% Q = B * uCov * B';
% R = mocapCov + D * uCov * D';  % trust on the equation of motion, given in Nm
% xk = T * [xv(1); (xv(2)-xv(1))*Fs];
% P = T * diag([mocapCov, 2*mocapCov*Fs^2]) * T';


Xs = smooth_LTI(A, B, C, D, Q, R, u, y, xk, P);

ry = Xs * C' + u * D' - y;
cost = mean(ry.^2) / R;  % cost function

