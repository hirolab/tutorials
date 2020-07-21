%% Generate measurements
clear
rng(10)

dlen = 1000 + 100;
Fs = 183;

stiff = 450;
damp = 2;
mass = 0.2;
sys = ss([0, 1; -stiff/mass, -damp/mass], [0; 1/mass], ...
         [1, 0; -stiff/mass, -damp/mass], [0; 1/mass]);  % KBM model, with acc measurement
     
time = (0 : dlen-1)' / Fs;
[bf, af] = fir1(10, 20/(Fs/2));
F = 10*filter(bf, af, randn(dlen, 1));  % band-limited white noise
[Y, ~, X] = lsim(sys, F, time);

forceCov = 0.5e0 ^ 2;  % force noise covariance
mocapCov = 2e-3^2;  % mocap noise covariance
noise = mvnrnd([0,0], [forceCov, 0; 0, mocapCov], dlen) * [1, 0*Fs^2; 0, 1]';  % some cross product op
Fv = F + noise(:,1);
xv = X(:,1) + noise(:,2);

offset = 10;
dlen = dlen - 100;

% admittance =====================
x0 = [xv(offset+1)
     (xv(offset+2) - xv(offset))*Fs/2];
P0 = diag([mocapCov, mocapCov*Fs]);
x = X(offset + (1:dlen)');
xd = X(offset + (1:dlen)',2);
xdd = Y(offset + (1:dlen)',2);
xv = xv(offset + (1:dlen)');
fv = Fv(offset + (1:dlen)');
f = F(offset + (1:dlen)');
time = time(offset + (1:dlen)');

figure, 
h1 = subplot(211); plot(time, [f, fv], '.-'), ylabel('Torque [Nm]')
h2 = subplot(212); plot(time, [x, xv] * 180/pi, '.-'), ylabel('Angle [deg]')
linkaxes([h1, h2], 'x')

% [xdf, xf, xddf] = deriv_sgolay(X(:,1), Fs, [3, 5]);
% figure, plot([X, Y(:,2)], '.-')
% hold on, set(gca, 'ColorOrderIndex', 1)
% plot([xf, xdf, xddf])

% figure, plot(xdd, '.-'), hold on, plot((1:dlen-1)'+0.5, diff(xd)*Fs, '.-')

%% Generate measurements (2 masses)
clear
rng(22)

dlen = 500 + 100;
Fs = 183;

stiff = 30;  
damp = 5;
mass = 0.02;
stiff2 = 500;  
damp2 = 5;
mass2 = 0.1;
mass3 = 1;
stiff3 = 100;
damp3 = 10;
sys = ss([0, 1, 0, 0, 0, 0
          [-stiff2, -damp2, stiff2, damp2, 0, 0]/mass2
          0, 0, 0, 1, 0, 0
          [stiff2, damp2, -(stiff2+stiff), -(damp2+damp), stiff, damp]/mass
          0, 0, 0, 0, 0, 1
          [0, 0, stiff, damp, -stiff-stiff3, -damp-damp3]/mass3], ...
          [0, 0; 1/mass2, 0; 0, 0; 0, -1/mass; 0, 0; 0, 1/mass3], ...
          [[-stiff2, -damp2, stiff2, damp2, 0, 0]/mass2
          [stiff2, damp2, -(stiff2+stiff), -(damp2+damp), stiff, damp]/mass
          [0, 0, stiff, damp, -stiff-stiff3, -damp-damp3]/mass3], ...
          [1/mass2, 0; 0, -1/mass; 0, 1/mass3]);  % KBM model, with acc measurement
     
time = (0 : dlen-1)' / Fs;
% [bf, af] = fir1(100, 20/(Fs/2));
[bf, af] = butter(5, 25/(Fs/2));
f0 = 0; f1 = 0.3/1;
F0 = polyval([f1, f0], time);
F = 4*filter(bf, af, randn(dlen, 1));  % band-limited white noise
[Y, ~, X] = lsim(sys, [F, F0], time);

forceCov = 5e-2 ^ 2;  % force noise covariance
mocapCov = 1e-3^2;  % mocap noise covariance
imuCov = 0.3e0 ^ 2;  % imu acc cov
noise = mvnrnd([0,0,0,0,0], diag([forceCov, mocapCov, mocapCov, mocapCov, imuCov]), dlen);  % some cross product op
Fv = F + noise(:,1);
xv = X(:,3) + noise(:,2);
x2v = X(:,1) + noise(:,3);
x3v = X(:,5) + noise(:,4);
x2ddv = Y(:,1) + noise(:,5);

offset = 10;
dlen = dlen - 100;

% admittance =====================
x = X(offset + (1:dlen)',3);
xd = X(offset + (1:dlen)',4);
xdd = Y(offset + (1:dlen)',2);
xv = xv(offset + (1:dlen)');
fv = Fv(offset + (1:dlen)');
f = F(offset + (1:dlen)');
x2 = X(offset + (1:dlen)',1);
x2d = X(offset + (1:dlen)',2);
x2dd = Y(offset + (1:dlen)',1);
x3 = X(offset + (1:dlen)',5);
x3d = X(offset + (1:dlen)',6);
x3dd = Y(offset + (1:dlen)',3);
x2v = x2v(offset + (1:dlen)');
x3v = x3v(offset + (1:dlen)');
time = time(offset + (1:dlen)');
x2ddv = x2ddv(offset + (1:dlen)');
F0 = F0(offset + (1:dlen)');

figure, 
h1 = subplot(311); plot(time, [f, fv], '.-'), ylabel('Torque [Nm]')
h2 = subplot(312); plot(time, [x, xv] * 180/pi, '.-'), ylabel('Angle [deg]')
h3 = subplot(313); plot(time, [x2dd, x2ddv] * 180/pi, '.-'), ylabel('Acc [deg/s^2]')
linkaxes([h1, h2, h3], 'x')

xn = [x2, x, x3];

% [xdf, xf, xddf] = deriv_sgolay(X(:,1), Fs, [3, 5]);
% figure, plot([X, Y(:,2)], '.-')
% hold on, set(gca, 'ColorOrderIndex', 1)
% plot([xf, xdf, xddf])

% rows = 10:dlen-10;
% xdz = (x(rows+1) - x(rows-1))*Fs/2;
% x3dz = (x3(rows+1) - x3(rows-1))*Fs/2;
% xddz = (x(rows+1) -2*x(rows) + x(rows-1))*Fs^2;
% x2ddz = (x2(rows+1) -2*x2(rows) + x2(rows-1))*Fs^2;

% this results in 0
% figure
% plot(mass*xdd + mass2*x2dd + damp*(xd-x3d) + stiff*(x-x3) - f + F0)

% figure
% plot([mass*xddz + mass2*x2ddz + damp*(xdz-x3dz) + stiff*(x(rows)-x3(rows)), f(rows)], '.-')


%% Solve with greybox
data = iddata(xv, fv, 1/Fs);
% data = iddata(zeros(dlen-2,1), [xv(3:end), fv(1:end-2)], 1/Fs);
% data = iddata(fv(1:end-2), [xv(3:end)], 1/Fs);

params = [500, 2, 0.1];  % [k, b]
aux = 1/Fs;  % [dt, m]

Minit = idgrey('build_spring_mass_damper', params, 'c', aux);

opt = greyestOptions;
opt.EnforceStability = true;
opt.DisturbanceModel = 'estimate';
opt.InitialState = 'estimate';
opt.Display = 'full';
Model = greyest(data,Minit,opt);

100*(Model.Report.Parameters.ParVector - [stiff; damp; mass]) ./ [stiff; damp; mass]
Model.Report.Parameters.ParVector'

figure, compare(data, Model, Minit)

%% LSQ

% select window and order depending on the signal bandwidth
[bf, af] = butter(5, 25/(Fs/2));
fvf = filtfilt(bf, af, fv);
xf = filtfilt(bf, af, xv);
[qd, qf, qdd] = deriv_sgolay(xf, Fs, [5, 7]);
x2f = filtfilt(bf, af, x2v);
[~, ~, q2dd] = deriv_sgolay(x2f, Fs, [5, 7]);
x3f = filtfilt(bf, af, x3v);
q3d = deriv_sgolay(x3f, Fs, [5, 7]);

% figure, plot(time, [qdd, xdd])

rows = length(bf) : (dlen-length(bf));  % remove filter transient
regressors = [xf-x3f, qd-q3d, qdd];
response = fvf - mass2*q2dd;

regressors(rows,1) = detrend(regressors(rows,1), 1);
response(rows,:) = detrend(response(rows,:), 1);

lm = fitlm(regressors(rows,:), response(rows,:));
figure, plot(lm)

params1 = [lm.Coefficients.Estimate(2:4); mass2]
100 * (lm.Coefficients.Estimate(2:4) - [stiff; damp; mass]) ./ [stiff; damp; mass]

%% Smart LSQ (removes correlation between regressors and samples)

[bf, af] = fir1(100, 25/(Fs/2));
xf = filtfilt(bf, af, xv);

rows = 20 : 20 : (dlen-20);  % remove filter transient
regressors = [xv(rows+7), xv(rows+6), xv(rows+5), xv(rows+4), xv(rows+3), xv(rows+2), xv(rows+1), xv(rows-0), ...
    xv(rows-1), xv(rows-2), xv(rows-3), xv(rows-4), xv(rows-5), xv(rows-6), xv(rows-7)];  % estimating impulse response
% regressors = [xf(rows+7), xf(rows+6), xf(rows+5), xf(rows+4), xf(rows+3), xf(rows+2), xf(rows+1), xf(rows-0), ...
%     xf(rows-1), xf(rows-2), xf(rows-3), xf(rows-4), xf(rows-5), xf(rows-6), xf(rows-7)];  % estimating impulse response
response = fv(rows);

lm = fitrlinear(regressors, response, 'OptimizeHyperparameters', 'all')
% xsol = lm.Coefficients.Estimate(2:end);
xsol = lm.Beta;
figure, plot((1:length(xsol))' - length(xsol)/2 - 0.5, xsol, 'o-'), xlabel('lag'), grid on

N = 10000;
xin = randn(N, 1);
yout = filter(xsol, 1, xin);
rows = 20 : (N-20);
xin = xin(rows);
yout = yout(rows + size(regressors,2)/2-0.5);
[Txy, freq] = tfestimate(xin, yout, hann(6*Fs), 3*Fs, 6*Fs, Fs);

s = 1j * 2*pi * freq;
Txy_non = mass*s.^2 + damp*s + stiff;

figure
h1 = subplot(211); semilogx(freq, 20*log10(abs([Txy, Txy_non])), 'linewidth', 1.5)
grid on, title('Magnitude'), ylabel('dB'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

h2 = subplot(212); semilogx(freq, angle([Txy, Txy_non])*180/pi, 'linewidth', 1.5)
grid on, title('Phase'), ylabel('deg'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

%% Frequency domain

[Txy, freq] = tfestimate(xv, fv, hann(2*Fs), 1*Fs, 2*Fs, Fs);
C = mscohere(xv, fv, hann(2*Fs), 1*Fs, 2*Fs, Fs);

s = 1j * 2*pi * freq;
Txy_non = mass*s.^2 + damp*s + stiff;

figure
h1 = subplot(311); semilogx(freq, 20*log10(abs([Txy, Txy_non])), '.-', 'linewidth', 1.5)
grid on, title('Magnitude'), ylabel('dB'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

h2 = subplot(312); semilogx(freq, angle([Txy, Txy_non])*180/pi, '.-', 'linewidth', 1.5)
grid on, title('Phase'), ylabel('deg'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

h3 = subplot(313); semilogx(freq, C, '.-', 'linewidth', 1.5)
grid on, title('Coherence'), ylabel('[1]'), xlabel('Frequency [Hz]')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

%% Transfer function

 % Import   data                       
data = iddata(xv, fv, 1/Fs);

% Transfer function estimation         
Options = tfestOptions;               
Options.Display = 'on';               
Options.WeightingFilter = [0, 30*2*pi];
                                       
 tf1 = tfest(data, 2, 0, Options);
 
fliplr(tf1.den / tf1.num)

%% moderate EM
rng(1)
% win=3, pol=2 ------------------------------------------------------------
% u = 0.5*(fv + fv([2:end, end]));  % shift force 0.5 sample
% uCov = forceCov * 2 * 0.5^2;
% % u = 0.5*(f + f([2:end, end])) + randn(dlen,1)*sqrt(forceCov);  % shift force 0.5 sample
% % uCov = forceCov;
% y = xv;
% T = diag([1e3, 1e3/Fs]);
% sys = c2d(ss([0, 1; -stiff/mass, -damp/mass], [0; 1/mass], [1, 0], 0), 1/Fs);
% A = T*sys.A/T; B = T*sys.B; C = sys.C/T; D = sys.D;
% Q = B * uCov * B';
% R = mocapCov + D * uCov * D';  % trust on the equation of motion, given in Nm
% xk = T * x0;
% P = T * P0 * T';

% forceCov = 0.5e0 ^ 2;  % force noise covariance
% mocapCov = 2e-3^2;  % mocap noise covariance

% 5th order ======================================
u = [xv(4:end), fv(1:end-3)];
uCov = diag([mocapCov, forceCov]);
y = zeros(size(u,1),1);
T = diag([1e3, 1e3, 1e3, 1e3, 1e3]); % eye(4); %diag([1, 1/Fs, 1/Fs^2]);  % state transformation
A = T * [0, 0, 0, 0, 0
         1, 0, 0, 0, 0
         0, 1, 0, 0, 0
         0, 0, 1, 0, 0
         0, 0, 0, 1, 0] / T;
B = T * [1, 0
         0, 0
         0, 0
         0, 0
         0, 0];
D = [0, 1];
Q = B * uCov * B'; %zeros(3);
R = 1e-3 ^ 2 + D * uCov * D';  % trust on the equation of motion, given in Nm
xk = T * [xv(1); xv(1); xv(1); xv(1); xv(1)];
P = T * diag(10*[mocapCov, mocapCov, mocapCov, mocapCov, mocapCov]) * T';
Cfun = @(params) -[params(1), params(2)*Fs, params(3)*Fs^2] * [[     0,   0,    1,    0,     0]
                                    [ -1/12, 2/3,    0, -2/3,  1/12]
                                    [ -1/12, 4/3, -5/2,  4/3, -1/12]] / T;

% % impedance D =======================
% u = [xv(3:end), fv(1:end-2)];
% uCov = diag([mocapCov, forceCov]);
% y = zeros(size(u,1),1);
% T = diag([1e3, 1e3, 1e3]); % eye(4); %diag([1, 1/Fs, 1/Fs^2]);  % state transformation
% A = T * [0, 0, 0
%          1, 0, 0
%          0, 1, 0] / T;
% B = T * [1, 0
%          0, 0
%          0, 0];
% D = [0, 1];
% Q = B * uCov * B'; %zeros(3);
% R = 1e-3 ^ 2 + D * uCov * D';  % trust on the equation of motion, given in Nm
% xk = T * [xv(1); xv(1); xv(1)];
% P = T * diag([mocapCov, mocapCov, mocapCov]) * T';
% 
% Cfun = @(params) -[params(1), params(2)*Fs, params(3)*Fs^2] * [[   0,  1,    0]
%                                     [ 1/2,  0, -1/2]
%                                     [   1, -2,    1]] / T;
                                
ny = size(y, 2);
N = size(y, 1);

params = [500; 0.5; 0.25];
lb = [1; 1e-3; 1e-3];
ub = [800; 10; 0.3];
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'UseParallel', false, 'MaxFunctionEvaluations', 3e4, ...
    'SpecifyObjectiveGradient', false, 'MaxIterations', 1e5);
    
Params = params;
Fval = NaN;
tic
for i = 1 : 400
    C = Cfun(params);
    [Xs, Ps, X, Ppost, E] = smooth_LTI(A, B, C, D, J, Q, R, u, y, xk, P);
    
    Sinv = zeros(ny, ny, N);  % error distro of residual
    for j = 1 : N
        Sinv(:,:,j) = inv(C * Ps(:,:,j) * C' + R);  % cov of e
    end
    resid_func = @(params) Xs * Cfun(params)' + u * D';
    fun = @(params) MLE_cost(params, resid_func, Sinv);
    
    [params, fval] = fmincon(fun, params, [], [], [], [], lb, ub, [], options);
    
    Params = [Params, params];
    Fval = [Fval, fval];
end
toc

figure
h1=subplot(411); plot(Params(1,:), '.-'), grid on, yline(stiff); ylim([0, 800]), title('Stiffness [Nm/rad]')
h2=subplot(412); plot(Params(2,:), '.-'), grid on, yline(damp); ylim([0, 10]), title('Damping [Nm/rad/s]')
h3=subplot(413); plot(Params(3,:), '.-'), grid on, yline(mass); ylim([0, 0.3]), title('Inertia [kg*m^2]')
h4=subplot(414); plot(Fval, '.-'), grid on, title('Cost')
linkaxes([h1, h2, h3, h4], 'x'), xlabel('iteration'), xlim([1, i])
set(gcf, 'position', [680   417   665   681])

params'
100 * (params - [stiff; damp; mass]) ./ [stiff; damp; mass]

%%

rows = (4:dlen-4)';
u = [x2v(rows + 3), xv(rows + 3), x3v(rows + 3)];
uCov = diag([mocapCov, mocapCov, mocapCov]);
y = [fv(rows), x2ddv(rows)];  % assume implicit function
yCov = diag([forceCov, imuCov]);

x0 = [x2v(rows(1)+2); x2v(rows(1)+1); x2v(rows(1)); x2v(rows(1)-1); x2v(rows(1)-2)
      xv(rows(1)+2); xv(rows(1)+1); xv(rows(1)); xv(rows(1)-1); xv(rows(1)-2)
      x3v(rows(1)+2); x3v(rows(1)+1); x3v(rows(1)); x3v(rows(1)-1); x3v(rows(1)-2)];
p0 = diag([mocapCov, mocapCov, mocapCov, mocapCov, mocapCov]);
P0 = blkdiag(p0, p0, p0);


obj = KBJ_impedance_implicit_3m(stiff, damp, mass, mass2, Fs, uCov, yCov);
[Xs, Ps, Xf, Pf] = smooth_NLTIimplicit(obj, u, y, x0, P0);


Xn = [x2(rows+2), x2(rows+1), x2(rows), x2(rows-1), x2(rows-2), ...
      x(rows+2), x(rows+1), x(rows), x(rows-1), x(rows-2), ...
      x3(rows+2), x3(rows+1), x3(rows), x3(rows-1), x3(rows-2)];
un = [x2(rows + 3), x(rows + 3), x3(rows + 3)];
yn = [f(rows), x2dd(rows)]; 


[goodnessOfFit(Xn + sqrt(mocapCov)*randn(size(Xn)), Xn, 'NRMSE'), ...
    goodnessOfFit(Xf, Xn, 'NRMSE'), goodnessOfFit(Xs, Xn, 'NRMSE')]

figure, plot([xv(rows,:), Xf(:,8), Xs(:,8), Xn(:,8)]), legend('measurement', 'filt', 'smooth', 'nominal')

%%

N = size(u, 1);
nx = length(x0);
nu = size(u, 2);
ny = size(y, 2);
nh = 2;

Xk = zeros(N, nx);
Hk = zeros(N, nh);

for i = 1 : N
    Xk(i,:) = StateTransition(obj, Xn(i,:)', un(i,:)');
    Hk(i,:) = Measurement(obj, Xn(i,:)', zeros(nu+ny,1), [un(i,:), yn(i,:)]');
end

goodnessOfFit(Xk(1:end-1,:), Xn(2:end,:), 'NRMSE')
figure, plot([yn+Hk, yn])


%% Test dynamic system object

% [obj, un, yn, x0n, P0n, Xn] = KBJ_impedance_implicit_3m_aug(stiff, damp, mass, mass2, Fs, mocapCov, forceCov, imuCov, ...
%         xn(:,1), xn(:,2), xn(:,3), f, x2dd);
% [obj, un, yn, x0n, P0n, Xn] = KBJ_impedance_implicit_3m_f0(stiff, damp, mass, mass2, ...
%     Fs, mocapCov, forceCov, x2, x, x3, f, 0, f0);
[obj, un, yn, x0n, P0n, Xn] = KBJ_impedance_implicit_3m_imu(stiff, damp, mass, mass2, ...
    Fs, 0, mocapCov, forceCov, imuCov, x2, x, x3, f, x2dd, F0);

Dlen = length(un);
X = zeros(Dlen, length(x0n));
Y = zeros(Dlen, 2);
for i = 1 : Dlen
    Y(i,:) = Measurement(obj, Xn(i,:)', zeros(3,1), [un(i,:), yn(i,:)]');
    X(i,:) = StateTransition(obj, Xn(i,:)', un(i,:)');
end 

% figure, plot([yn+Y, yn], '.-')
% figure, plot([yn+Y-yn], '.-'), ylim([-2, 2])
% figure, plot([X(1:end-1,[3, 8, 13]), Xn(2:end,[3, 8, 13])], '.-')
% figure, plot([X(1:end-1,[3, 8, 13])-Xn(2:end,[3, 8, 13])], '.-')

%%

rows = (10:dlen-10)';

% params = [300; 1; mass; mass2];
% params = params1;
params = [stiff; damp; mass; mass2];
[obj, u, y, x0, P0, x] = KBJ_impedance_implicit_3m_aug(params(1), params(2), params(3), params(4), Fs, mocapCov, forceCov, imuCov, ...
        x2v, xv, x3v, fv, x2ddv);
% [obj, u, y, x0, P0, x] = KBJ_impedance_implicit_3m_aug(params(1), params(2), params(3), params(4), Fs, mocapCov, forceCov, imuCov, ...
%         xn(:,1), xn(:,2), xn(:,3), f, x2dd);
[Xs, Ps, Xf, Pf] = smooth_NLTIimplicit(obj, u, y, x0, P0);

Dlen = length(rows);
Yf = zeros(Dlen, 2);
Ys = zeros(Dlen, 2);
for i = 1 : Dlen
    Yf(i,:) = Measurement(obj, Xf(i,:)', zeros(2,1), [u(i,:), y(i,:)]');
    Ys(i,:) = Measurement(obj, Xs(i,:)', zeros(2,1), [u(i,:), y(i,:)]');
end

% figure, plot([y, y+Yf, y+Ys, yn])
figure, plot([y, yn], '.-'), grid on
figure, plot([y+Yf, yn], '.-'), grid on


Xs = Xs / obj.T';
Xf = Xf / obj.T';

% [goodnessOfFit([x2v(rows), xv(rows), x3v(rows)], xn(rows,:), 'NRMSE'), ...
%     goodnessOfFit(Xf(:,[3,8,12]),xn(rows,:), 'NRMSE'), ...
%     goodnessOfFit(Xs(:,[3,8,12]), xn(rows,:), 'NRMSE')]

%     figure, plot([x2v(rows), Xs(:,3), xn(rows,1)])
%     figure, plot([xv(rows), Xs(:,8), xn(rows,2)])
%     figure, plot([x3v(rows), Xs(:,13), xn(rows,3)])



Params = Xs(:,end-3:end);
Paramsf = Xf(:,end-3:end);
figure, set(gcf, 'position', [680          94         637        1004])
h1=subplot(411); plot(time(rows), [Params(:,1), Paramsf(:,1)], '.-'), grid on, yline(stiff); ylim([0, 2*stiff]), title('Stiffness [Nm/rad]')
h2=subplot(412); plot(time(rows), [Params(:,2), Paramsf(:,2)], '.-'), grid on, yline(damp); ylim([0, 2*damp]), title('Damping [Nm/rad/s]')
h3=subplot(413); plot(time(rows), [Params(:,3), Paramsf(:,3)], '.-'), grid on, yline(mass); ylim([0, 2*mass]), title('Inertia [kg*m^2]')
h4=subplot(414); plot(time(rows), [Params(:,4), Paramsf(:,4)], '.-'), grid on, yline(mass2); ylim([0, 2*mass2]), title('Inertia 2 [kg*m^2]')
% h5=subplot(515); plot(Fval, '.-'), grid on, title('Cost')
linkaxes([h1, h2, h3, h4], 'x'), xlabel('iteration'),% xlim([1, length(Fval)])
drawnow()

params = mean(Xs(:,end-3:end))';

params'
100 * (params - [stiff; damp; mass; mass2]) ./ [stiff; damp; mass; mass2]


%% moderate EM (rates are states)
% rng(1)

% params = [24; 2.3; 0.025; mass2];
params = params1;
% params = [stiff; damp; mass; mass2];
lb = [1; 1e-3; 1e-2; mass2];
ub = [800; 10; 0.3; mass2];
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'UseParallel', false, 'MaxFunctionEvaluations', 3e4, ...
    'SpecifyObjectiveGradient', true, 'MaxIterations', 1e5, 'CheckGradients', false);

[~, un, yn, ~, ~, Xn] = KBJ_impedance_implicit_3m_imu(0,0,0,0, Fs, 0, mocapCov, forceCov, imuCov, x2, x, x3, f, x2dd, F0);
% [~, un, yn, ~, ~, Xn] = KBJ_impedance_implicit_3m_f0(0,0,0,0, Fs, mocapCov, forceCov, x2, x, x3, f, 0, f0);
% [~, un, yn, ~, ~, Xn] = KBJ_impedance_implicit_3m(0,0,0,0, Fs, mocapCov, forceCov, x2, x, x3, f, 0);

f1 = figure(); set(gcf, 'position', [680          76         665        1022])
Params = params;
Fval = NaN;
tic
for i = 2 : 10000
    % params = [stiff; damp; mass; mass2];
    resCov = 0.1*exp(-i/50);

    [obj, u, y, x0, P0, X] = KBJ_impedance_implicit_3m_imu(params(1), params(2), params(3), params(4), ...
        Fs, resCov, mocapCov, forceCov, imuCov, x2v, xv, x3v, fv, x2ddv, xv*0);
%     [obj, u, y, x0, P0, X] = KBJ_impedance_implicit_3m_f0(params(1), params(2), params(3), params(4), ...
%         Fs, mocapCov, forceCov, x2v, xv, x3v, fv, resCov, 0);
%     [obj, u, y, x0, P0, X] = KBJ_impedance_implicit_3m(params(1), params(2), params(3), params(4), ...
%         Fs, mocapCov, forceCov, x2v, xv, x3v, fv, resCov);
    [Xs, Ps, Xf, Pf] = smooth_NLTIimplicit(obj, u, y, x0, P0);

%     figure, plot([x2v(rows), Xs(:,3), Xn(:,3)])
%     figure, plot([xv(rows), Xs(:,8), Xn(:,8)])
%     figure, plot([x3v(rows), Xs(:,13), Xn(:,13)])
%     goodnessOfFit(Xs(:,[3,8,12]), Xn(:,[3,8,12]), 'NRMSE')
%     goodnessOfFit(Xf(:,[3,8,12]), Xn(:,[3,8,12]), 'NRMSE')
%     goodnessOfFit([x2v(rows), xv(rows), x3v(rows)], [x2(rows), x(rows), x3(rows)], 'NRMSE')
%     figure, plot([squeeze(Ppost(8,8,:)), squeeze(Ps(8,8,:))], '.-')
% 
%     figure, plot([xv(rows), Xf(:,8), Xs(:,8), Xn(:,8)]), legend('noise', 'filt', 'smooth', 'nominal')

    [Sinv, Q, func] = wrap_cost_fcn(obj, Xs, Ps, u, y);
    fun = @(params) MLE_cost(params, func, Sinv, Q);
    [params, fval, exitflag, output] = fmincon(fun, params, [], [], [], [], lb, ub, [], options);
    
%     opts = optimoptions(@fmincon, ...
%         'Algorithm', 'interior-point', ...
%         'StepTolerance', 1e-10, ...
%         'MaxFunctionEvaluations', 3e4, ...
%         'Display','off', ...
%         'SpecifyObjectiveGradient', false);
%     problem = createOptimProblem('fmincon', ...
%         'objective', fun, ...
%         'x0', params, ...
%         'lb', lb, ...
%         'ub', ub, ...
%         'Aineq', [], ...
%         'bineq', [], ...
%         'options',opts);
%     ms = MultiStart('Display', 'off', 'UseParallel', false);
%     tic
%     [params, fval] = run(ms, problem, 30);
%     toc
    
    Params = [Params, params];
    Fval = [Fval, fval];
    
    if mod(i, 10) == 0 || all(abs(Params(:,end-1) - params) ./ params == 0)
        h1=subplot(511); plot(Params(1,:), '.-'), grid on, yline(stiff); ylim([0, max(2*stiff, max(Params(1,:)))]), title('Stiffness [Nm/rad]')
        h2=subplot(512); plot(Params(2,:), '.-'), grid on, yline(damp); ylim([0, max(2*damp, max(Params(2,:)))]), title('Damping [Nm/rad/s]')
        h3=subplot(513); plot(Params(3,:), '.-'), grid on, yline(mass); ylim([0, max(2*mass, max(Params(3,:)))]), title('Inertia [kg*m^2]')
        h4=subplot(514); plot(Params(4,:), '.-'), grid on, yline(mass2); ylim([0, 2*mass2]), title('Inertia 2 [kg*m^2]')
        h5=subplot(515); plot(Fval, '.-'), grid on, title('Cost')
        linkaxes([h1, h2, h3, h4, h5], 'x'), xlabel('iteration'), xlim([1, length(Fval)])
        drawnow()
    end
    
    if all(abs(Params(:,end-1) - params) ./ params == 0) && i > 20
        break
    end
end
toc

% ord = 5;
% figure, plot([xv(rows), Xs(:,ord+(ord+1)/2), Xn(:,ord+(ord+1)/2)], '.-')
% figure, plot([x2v(rows), Xs(:,(ord+1)/2), Xn(:,(ord+1)/2)], '.-')

params'
100 * (params - [stiff; damp; mass; mass2]) ./ [stiff; damp; mass; mass2]


%% Kalman Smoother
rng(1)
% win=3, pol=2 ------------------------------------------------------------
% u = 0.5*(fv + fv([2:end, end]));  % shift force 0.5 sample
% uCov = forceCov * 2 * 0.5^2;
% % u = 0.5*(f + f([2:end, end])) + randn(dlen,1)*sqrt(forceCov);  % shift force 0.5 sample
% % uCov = forceCov;
% y = xv;
% T = diag([1e3, 1e3/Fs]);
% sys = c2d(ss([0, 1; -stiff/mass, -damp/mass], [0; 1/mass], [1, 0], 0), 1/Fs);
% A = T*sys.A/T; B = T*sys.B; C = sys.C/T; D = sys.D;
% Q = B * uCov * B';
% R = mocapCov + D * uCov * D';  % trust on the equation of motion, given in Nm
% xk = T * x0;
% P = T * P0 * T';

% % impedance D =======================
% u = [xv(3:end), fv(1:end-2)];
% uCov = diag([mocapCov, forceCov]);
% y = zeros(size(u,1),1);
% nbuffer = 3;
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
% xk = T * [xv(1); xv(1); xv(1)];
% P = T * diag([mocapCov, mocapCov, mocapCov]) * T';

% 5th order
u = [xv(4:end), fv(1:end-3)];
uCov = diag([mocapCov, forceCov]);
y = zeros(size(u,1),1);
T = diag([1e3, 1e3, 1e3, 1e3, 1e3]); % eye(4); %diag([1, 1/Fs, 1/Fs^2]);  % state transformation
A = T * [0, 0, 0, 0, 0
         1, 0, 0, 0, 0
         0, 1, 0, 0, 0
         0, 0, 1, 0, 0
         0, 0, 0, 1, 0] / T;
B = T * [1, 0
         0, 0
         0, 0
         0, 0
         0, 0];
C = -[stiff, damp*Fs, mass*Fs^2] * [[     0,   0,    1,    0,     0]
                                    [ -1/12, 2/3,    0, -2/3,  1/12]
                                    [ -1/12, 4/3, -5/2,  4/3, -1/12]] / T;
D = [0, 1];
Q = B * uCov * B'; %zeros(3);
R = 1e-3 ^ 2 + D * uCov * D';  % trust on the equation of motion, given in Nm
xk = T * [xv(1); xv(1); xv(1); xv(1); xv(1)];
P = T * diag(10*[mocapCov, mocapCov, mocapCov, mocapCov, mocapCov]) * T';



rows = (10:dlen-10)';
xx = [x(rows+2), x(rows+1), x(rows-0), x(rows-1), x(rows-2)] * T';
% xx = [x(rows+1), x(rows-0), x(rows-1)] * T';
figure, plot([-xx * C'- f(rows)], '.-'), grid on, ylim([-0.2, 0.2])

return

[Xs, Ps, X, Ppost, E] = smooth_LTI(A, B, C, D, Jy, Q, R, u, y, xk, P);



% figure, plot(x(3:end), Xs * C' + u * D', '.')
ry = Xs * C' + u * D' - y;
sum(R * ry.^2)  % cost function
% 
% return

X = X / T';
Xs = Xs / T';

% figure, plot([X(:,1), Xs(:,1)], '.-')
% hold on, plot((1:dlen)'+-1, x(1:end), '--')
% plot((1:dlen)'+-1, xv(1:end), 'o-'), 
% legend('filtered', 'smoothed', 'nominal', 'measurement')
% 
% figure, plot([X(:,1), Xs(:,1)] - x, '.-')
% legend('filtered', 'smoothed'), grid on

% mean([X(:,1), Xs(:,1)] - x(2:end-1))
% std([X(:,1), Xs(:,1)] - x(2:end-1))

% mean([X(:,1), Xs(:,1)] - x(1:end-1))
% std([X(:,1), Xs(:,1)] - x(1:end-1))

% mean([X(:,1), Xs(:,1)] - x(1:end))
% std([X(:,1), Xs(:,1)] - x(1:end))
% figure, autocorr(Xs(:,1) - x(1:end))


% ===========
% figure, plot([Xs(:,2),x(1:end-2)], '.-'), ylim([-1e-3, 1e-3])  % ord 3
% figure(1), hold on, plot([Xs(:,3)-x(1:end-3)], '.-'), ylim([-1e-3, 1e-3])  % ord 5
%=========

% hold on, plot((1:dlen)'+0, f, '--')
% legend('smoothed', 'nominal')

% 
% figure, plot([X(:,end), Xs(:,end)] - f, '.-')
% legend('filtered', 'smoothed'), grid on
% 
% figure, autocorr(E)
% fprintf('%.3f %c %.3f N\n', mean(E), 177, std(E))

%%
rows= 10:dlen-10;

% figure, plot([(x(rows+1)-x(rows-1))*Fs/2, xd(rows)], '.-')
% figure, plot([(x(rows+1)-2*x(rows)+x(rows-1))*Fs^2, xdd(rows)], '.-')

[xd2, xf2, xdd2] = deriv_sgolay(x, Fs, [2, 3]);
% figure, plot([xd2(rows), xd(rows)], '.-')
figure, plot([xdd2(rows), xdd(rows)], '.-')




%%

nx = 4;
nu = 2;
ny = 1;

psi = zeros(nx+ny, nx+nu);
sig = zeros(nx+nu, nx+nu);
for i = 1 : dlen-1
    psi = psi + [Xs(i+1,:)' * Xs(i,:) + M(:,:,i+1), Xs(i+1,:)' * u(k,:)
                 z(i,:)' * Xs(i,:), z(i,:)' * u(k,:)] / (dlen-1);
    sig = sig + [Xs(i,:)' * Xs(i,:) + Ps(:,:,i), Xs(i,:)' * u(i,:)
                 u(i,:)' * Xs(i,:), u(i,:)' * u(i,:)] / (dlen-1);
end

sol = psi / sig;
sol

AA = sol(1:nx, 1:nx);
BB = sol(1:nx, nx+1:end);
CC = sol(nx+1:end, 1:nx);
DD = sol(nx+1:end, nx+1:end);

%% Generate C function (time sequence to q, qd, qdd)

% given N points, what is x, xd, and xdd at center?
q = sym('q', [7, 1]);
Fs = 1;
time = (3:-1:-3)' / Fs;
syms q0 qd qd2 qd3 qd4 qd5 qd6
% eqs = q == q0 + time*qd + qd2*time.^2/2 + qd3*time.^3/6 + qd4*time.^4/24 + qd5*time.^5/120 + qd4*time.^6/720;

qderiv = linsolve([ones(7,1), time, time.^2/2, time.^3/6, time.^4/24, time.^5/120, time.^6/720], q);
[A,b] = equationsToMatrix(qderiv(1:3), q)

%%
n = 11;
time = ((n-1)/2 : -1 : -(n-1)/2)' / Fs;
% (time.^(0:(n-1)) ./ factorial(0:(n-1))) 
deriv_blk = inv((time.^(0:(n-1)) ./ factorial(0:(n-1))))






