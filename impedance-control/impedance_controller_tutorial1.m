%% Define sim params
dt = 0.05;
time = (0 : dt : 5)';  % simulation interval
fext_fcn = @(t) 10*(t < 1);  % external disturbance
X0 = [0; 0];  % initial condition (displacement, speed)

figure, plot(time, fext_fcn(time), '.-')
ylabel('external force [N]'), xlabel('time [s]')

%% Original system
m = 1;
b = 1;
k = 10;

u_fcn = @(t) zeros(size(t));
f = @(t, x) [0, 1; -k/m, -b/m] * x + [0; 1/m] * (fext_fcn(t) + u_fcn(t));
[~, x] = ode45(f, time, X0);
y = x * [1, 0]';

figure, plot(time, y, '.-')
ylabel('displacement [m]'), xlabel('time [s]')
title('Original')

%% Sim desired system
%   ignore desired trajectory for now

ms = 0.1;
bs = 0.1;
ks = 10;

u_fcn = @(t) zeros(size(t));
f = @(t, x) [0, 1; -ks/ms, -bs/ms] * x + [0; 1/ms] * (fext_fcn(t) + u_fcn(t));
[~, x] = ode45(f, time, X0);
y = x * [1, 0]';

figure, plot(time, y, '.-')
ylabel('displacement [m]'), xlabel('time [s]')
title('Desired')

%% Sim impedance controller (closed loop system)

u_fcn = @(t, x, fext) (m/ms)*(-[ks, bs]*x + fext) + [k, b]*x - fext;  % ctrl law
f = @(t, x) [0, 1; -k/m, -b/m] * x + [0; 1/m] * (fext_fcn(t) + u_fcn(t, x, fext_fcn(t)));
[~, x] = ode45(f, time, X0);
y = x * [1, 0]';

hold on, plot(time, y, '.-')
ylabel('displacement [m]'), xlabel('time [s]')
title('Closed-Loop')

