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


%%

