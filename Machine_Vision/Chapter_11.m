cd('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training')
clc; clear all; close all;

%% remember to get the toolkits from:
% https://petercorke.com/toolboxes/robotics-toolbox/
% https://petercorke.com/toolboxes/machine-vision-toolbox/

%% Book title: Robotics, Vision and Control Fundamental algorithms in MATLAB: Second Edition pdf
% Can be accessed via Purdue Library

%% 11.1.2 Modeling a perspective Camera
%Principal point is synonymous to image center, where the 2D image is
%projected

cam = CentralCamera('focal', 0.015); %returns an instance of a CentralCamera object with a 15 mm lens

% if we define the following world point.
P = [0.3, 0.4, 3.0]';
cam.project(P) %output is int on image plane

cam.project(P, 'pose', SE3(-0.5, 0, 0) )% moving camera to the left

%% 11.1.3 Discrete Image Plane

cam = CentralCamera('focal', 0.015, 'pixel', 10e-6,	'resolution', [1280 1024], 'centre', [640 512], 'name', 'mycamera')
P = [0.3, 0.4, 3.0]';
cam.project(P)


%% 11.1.4 Camera Matrix
%The camera intrinsic parameter matrix K for this camera is:
cam.K

%camera matrix is implicitly created when the Toolbox camera object is constructed
%and for this example is
cam.C

%camera field of view  
cam.fov() * 180/pi  % in degrees

%% 11.1.5 Projecting Points
 
P = mkgrid(3, 0.2, 'pose', SE3(0, 0, 1.0)); %needs a correction in one of the functions to work
%modify ES3.convert to ES3.check

P(:,1:4)
cam.project(P)
cam.plot(P)
%%
Tcam = SE3(-1,0,0.5)*SE3.Ry(0.9);

%results in an oblique view of the plane
cam.plot(P, 'pose', Tcam)

cam.project([1 0 0 0]', 'pose', Tcam);
p = cam.plot(P, 'pose', Tcam) %this code doesn't work :(

%%
cube = mkcube(0.2, 'pose', SE3(0, 0, 1) );
cam.plot(cube);
[X,Y,Z] = mkcube(0.2, 'pose', SE3(0, 0, 1), 'edge');
cam.mesh(X, Y, Z)


Tcam = SE3(-1,0,0.5)*SE3.Ry(0.8);
cam.mesh(X, Y, Z, 'pose', Tcam);
%% This one did work! And  very cool plot

theta = [0:500]/100*2*pi;
[X,Y,Z] = mkcube(0.2, [], 'edges');
for th=theta
T_cube = SE3(0, 0, 1.5)*SE3.rpy(th*[1.1 1.2 1.3])
cam.mesh( X, Y, Z, 'objpose', T_cube ); drawnow
end


%% 11.1.6 Lens Distortion

cam = CentralCamera('focal', 0.015, 'pixel', 10e-6, ...
'resolution', [1280 1024], 'centre', [512 512], ...
'distortion', [k1 k2 k3 p1 p2] )


%% 11.2.1 Camera Calibration
P = mkcube(0.2);
T_unknown = SE3(0.1, 0.2, 1.5) * SE3.rpy(0.1, 0.2, 0.3);
cam = CentralCamera('focal', 0.015, ...
'pixel', 10e-6, 'resolution', [1280 1024], 'noise', 0.05);
p = cam.project(P, 'objpose', T_unknown);

C = camcald(P, p)

%maximum residual is expressed and matric C is the camera matrix
%% Chapter 11.2.2  Decomposing the Camera Calibration Matrix
null(C)' % is the world origin of the camera frame
h2e(ans)' %expressed in cartesian form 
T_unknown.inv.t' %close to the true value
est = invcamcal(C) %recover other camera parameters from camera matrix WOW!
%which returns a CentralCamera object with its parameters set to values that result
%in the same camera matrix. We note immediately that the focal length is very large
%compared to the true focal length of our lens which was 0.015 m, and that the pixel
% sizes are very large. From Eq. 11.9 we see that focal length and pixel dimensions always
% appear together as factors f /?w and f /? h. The function invcamcal has set
% ?w= 1 but the ratios of the estimated parameters

est.f/est.rho(1) % ratio comparison of estimate
cam.f/cam.rho(2) % and true

% this following code shows the error between the estimates in the camera
% world position
hold on; plot_sphere(P, 0.03, 'r')
trplot(eye(4,4), 'frame', 'T', 'color', 'b', 'length', 0.3)

est.plot_camera()

%% 11.2.3 Pose Estimation

cam = CentralCamera('focal', 0.015, 'pixel', 10e-6, ...
'resolution', [1280 1024], 'centre', [640 512]); % camera with known paramers

P = mkcube(0.2); % we want to determine the pose of the cube

T_unknown = SE3(0.1, 0.2, 1.5) * SE3.rpy(0.1, 0.2, 0.3);
T_unknown.print
p = cam.project(P, 'objpose', T_unknown);

T_est = cam.estpose(P, p).print

%% 11.2.4 Camera Calibration Toolbox
addpath('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training\TOOLBOX_calib')
calib_gui  %Launches Calibration Toolbox GUI

















































