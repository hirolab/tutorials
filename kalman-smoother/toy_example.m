%% Toy example first
%% Create data
clear, clc

dt = 1/350;
t = (0:dt:0.5)';
u = zeros(length(t), 0); %mvnrnd(0, dt*1, length(t));
x0 = [0.01; 0.9];
A = expm([0, 1; -100/0.1, -3/0.1]*dt);
% C = [1, 0; 0.5, 0; -1.2, 0];
C = [1, 0];
sysd = ss(A, zeros(2,0), C, zeros(1, 0), dt);

[y, ~, x] = lsim(sysd, u, t, x0);
R = dt*diag([1e-4, 0.9e-4, 1.1e-4]);
z = y + mvnrnd(0, R(1,1), length(t));
% figure, plot([y, z])

figure
plot(t, [y, z])

%% Constant acceleration model

% A = expm([0, 1, 0; 0, 0, 1; 0, 0, 0]*dt);
% C = [1, 0, 0; 0.5, 0, 0; -1.2, 0, 0];

A = expm([0, 1; -100/0.1, -3/0.1]*dt);
C = [1, 0];

z = z(end:-1:1, :);
x = x(end:-1:1,:) .* [1, -1];

data = iddata(z, u, dt);

[xdh1, xh1] = deriv_sgolay(z, 1/dt, [5, 21]);  % traditional method
xh1 = [xh1(:,1), xdh1(:,1)];


%% Filter using the correct model
% Cons: Due to the wrong model assumption of constant acceleration, the
% filter has to update the acceleration and lags behind. A Smoother should
% fix this problem

Q = [   dt^5/20,   dt^4/8,  dt^3/6
        dt^4/8,     dt^3/3, dt^2/2
        dt^3/6,     dt^2/2, dt];  % white noise in accel
initialStateGuess = [z(1,1)
                    (z(2,1) - z(1,1)) / dt
                    (z(3,1) - 2*z(2,1) + z(1,1)) / dt^2];
    
    % Q = [   dt^3/3, dt^2/2
%         dt^2/2, dt];  % white noise in vel
% initialStateGuess = [z(1,1); (z(2,1)-z(1,1))/dt];

state_fcn = @(X) A*X;
measurement_fcn = @(X) C*X;

kf = extendedKalmanFilter(state_fcn, ...
                          measurement_fcn, ...
                          initialStateGuess, ...
                          'ProcessNoise', Q * 1e2, ...
                          'MeasurementNoise', R, ...
                          'StateCovariance', Q * 1e5);

num_samples = length(t);
xh = zeros(num_samples, length(initialStateGuess));
yh = zeros(num_samples, size(y,2));
for k = 1 : num_samples
    [CorrectedState,CorrectedStateCovariance] = correct(kf, z(k,:));
    xh(k,:) = CorrectedState;
    yh(k,:) = measurement_fcn(CorrectedState);
    
    [PredictedState,PredictedStateCovariance] = predict(kf);
end

figure, 
h1 = subplot(221); plot(t, x(:,1), '--', t, xh(:,1), '.-', t, xh1(:,1), '.'), grid on
h2 = subplot(223); plot(t, x(:,2), '--', t, xh(:,2), '.-', t, xh1(:,2), '.'), grid on
legend('reference', 'smoothed', 'numerical')
linkaxes([h1, h2], 'x')

rows = (t > 0.1) & (t < (t(end) - 0.1));
subplot(222), autocorr(x(rows,1) - xh(rows,1))
subplot(224), autocorr(x(rows,2) - xh(rows,2))
sgtitle('Kalman Filter')

%% Smoothing via greybox

par1 = 1e0;  % Initial guess for x0_0
par2 = z(1);  % Initial guess for x0_0
par3 = (z(2) - z(1))/dt;
% par4 = (z(3) - 2*z(2) + z(1)) / dt^2;
Pvec = [par1; par2; par3];
auxVal = [1e5, dt]; % R2=1

Minit = idgrey('mynoise',Pvec,'d',auxVal);
Minit.Structure.Parameters.Free(1) = false;

opt = greyestOptions;
opt.EnforceStability = true;
opt.DisturbanceModel = 'model';
% opt.Focus = 'simulation';
opt.InitialState = 'estimate';
opt.Display = 'full';
Model = greyest(data,Minit,opt);

% figure, compare(data, Model)
yp = predict(Model, data, 1);
yh = yp.OutputData(:,1);

xdh = deriv_sgolay(yh, 1/dt, [5, 21]);
xh = [yh, xdh];

figure, 
h1 = subplot(221); plot(t, x(:,1), '--', t, xh(:,1), '.-', t, xh1(:,1), '.'), grid on
h2 = subplot(223); plot(t, x(:,2), '--', t, xh(:,2), '.-', t, xh1(:,2), '.'), grid on
legend('reference', 'smoothed', 'numerical')
linkaxes([h1, h2], 'x')

rows = (t > 0.1) & (t < (t(end) - 0.1));
subplot(222), autocorr(x(rows,1) - xh(rows,1))
subplot(224), autocorr(x(rows,2) - xh(rows,2))
sgtitle('Grey Box')

%% Predict 

[yh, ~, xh] = lsim(Model, u, t, Model.Report.Parameters.X0);
% xdh = deriv_sgolay(yh, 1/dt, [5, 11]);
xdh = deriv_sgolay(z, 1/dt, [5, 15]);
figure, 
h1 = subplot(211); plot(t, x(:,1), '--', t, xh(:,1), '-'), grid on
h2 = subplot(212); plot(t, x(:,2), '--', t, xh(:,2), '-', t, xdh, '.-'), grid on
linkaxes([h1, h2], 'x')

%% Smooth

xh = zeros(num_samples, 3);
yh = zeros(num_samples, 1);

A = [   1,  dt, dt^2
        0,  1,  dt
        0,  0,  1];  % constant acceleration
C = [1, 0, 0];
initialStateGuess = [z(1); (z(2)-z(1))/dt; (z(3)-2*z(2)+z(1))/(dt^2)];
sysEst = ssm(A, zeros(3,0), C, 'Mean0', initialStateGuess,'StateType',[2;2;2]);

xh = smooth(sysEst, z);

figure, plot(t, xh(:,1:2), '.-')
hold on, set(gca, 'ColorOrderIndex', 1)
plot(t, x)

%% Smooth latent state







