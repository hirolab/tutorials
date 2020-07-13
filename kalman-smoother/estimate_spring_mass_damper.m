%% Estimate MRF impedance
%% Create artificial data
clear, clc

dt = 1/180;
t = (0:dt:1.5)';
u = zeros(length(t), 0); % no input
x0 = [2e-4; 9e-3];
A = ([0, 1; -4000/1.1, -5/1.1]);
C = [1, 0];
sys = ss(A, zeros(2,0), C, zeros(1, 0));

[y, ~, x] = lsim(sys, u, t, x0);
R = dt*3e-8;
z = y + mvnrnd(R(1,:)*0, R, length(t));

data = iddata(z, u, dt);

figure
plot(t, [y, z])

%% Load experimental data
clear, clc
qtmp = QTMParser('C:/Users/garamizo/Downloads/MRF_Test1_NoMagnets.mat');

dt = 1 / qtmp.Fs;

figure
% fplate
subplot(231)
QTMParser.plotBodyMarkers(qtmp.bodies(1), qtmp.markers(1:4))

% pin-holder
subplot(232)
QTMParser.plotBodyMarkers(qtmp.bodies(2), qtmp.markers(5:8))

% beam
subplot(233)
QTMParser.plotBodyMarkers(qtmp.bodies(3), qtmp.markers(9:12))

% L-frame
subplot(234)
QTMParser.plotBodyMarkers(qtmp.bodies(5), qtmp.markers(13:16))

% weight2
subplot(235)
QTMParser.plotBodyMarkers(qtmp.bodies(6), qtmp.markers(17:19))

% select "Beam - 3" marker z axis as input
rows = 7437 + (1:300);
z = detrend(qtmp.markers(11).trans(rows,3), 0);
figure, plot(z), grid on

data = iddata(z, [], dt);

%% Smoothing via greybox

params = [4000, 3];  % [k, b]
aux = [dt, 1.1];  % [dt, m]

Minit = idgrey('build_spring_mass_damper', params, 'd', aux);

opt = greyestOptions;
opt.EnforceStability = true;
opt.DisturbanceModel = 'estimate';
opt.InitialState = 'estimate';
opt.Display = 'full';
Model = greyest(data,Minit,opt);

figure, compare(data, Model)

