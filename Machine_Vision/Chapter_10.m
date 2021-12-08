%% Chapter 10 Light and Color

clc; clear all; close all;

%% remember to get the toolkits from:
% https://petercorke.com/toolboxes/robotics-toolbox/
% https://petercorke.com/toolboxes/machine-vision-toolbox/

%% Book title: Robotics, Vision and Control Fundamental algorithms in MATLAB: Second Edition pdf
% Can be accessed via Purdue Library

%%
lambda = [300:10:1000]*1e-9
%lambda = 500:100:600 %wavelength 
sun = blackbody(lambda, 5778);

for T=3000:3000:6000
plot( lambda, sun); hold all
legend('temperature 3000','6000')
end
%%
human = luminos(lambda);
plot(lambda, human)

%% 10.2.3 Reproducing Colors
%lambda is the color wavelength

green = cmfrgb(500e-9) % this returns how much RGB is needed to return a 500 color wavelength (green)
%however, we cannot have negative red
% so we add some white to make it possible
white = -min(green) * [1 1 1]
feasible_green = green + white
%% 10.2.4 Chromaticity Space

%The tristimulus values describe color
% r = R/(R+G+B) and so on for blue and gree

%lambda is the color wavelength

% the following plots the color of the rainbow (human visual ranges)
[r,g] = lambda2rg( [400:700]*1e-9 );
figure;plot(r, g)
rg_addticks

primaries = lambda2rg( cie_primaries() );
figure;plot(primaries(:,1), primaries(:,2), 'o')

% The COLOR GAMUT is the space where colors can be created. Colors made
% from primaries can only be created within the ranges of those colors.
% That is why low wavelength light cannot be created by primaries.
% another limitation is that some colors require negative amounts of a
% color however this is not possible and we can only add so much light
% before being limited.

green_cc = lambda2rg(500e-9)
figure;plot2(green_cc, 's')
white_cc = tristim2cc([1 1 1])
plot2(white_cc, 'o')

cmf = cmfxyz(lambda);
figure;plot(lambda, cmf);

% CIE XYZ primary that is required to match a spectral color and we note that these curves are never negative.
[x,y] = lambda2xy(lambda);
figure;plot(x, y);
showcolorspace('xy')

%% 10.2.6 Other Color and Chromaticity Spaces
colorspace('RGB->HSV', [1, 0, 0])
colorspace('RGB->HSV', [0, 1, 0])
colorspace('RGB->HSV', [0, 0, 1])

flowers = iread('flowers4.png', 'double');
hsv = colorspace('RGB->HSV', flowers); %changed color space from RGB to HSV
about hsv

colorspace('RGB->HSV', [0, 0.5, 0] + [0.4, 0.4, 0.4]) %green hue and medium saturatio value

figure;idisp( hsv(:,:,1) )
figure;idisp( hsv(:,:,2) )
figure;idisp( hsv(:,:,3) )

Lab = colorspace('RGB->Lab', flowers);
about Lab

idisp(Lab)
idisp( Lab(:,:,2) )
idisp( Lab(:,:,3) )

%% 10.3.2

% For example at night under a yellowish tungsten lamp the pages of a book still
% appear white to us, but a photograph of that scene viewed later under different lighting
% conditions will look yellow. All of this poses real problems for a robot that is using color to
% understand the scene because the observed chromaticity varies with lighting

% the following segment compares how a red brick can look different under
% two different light sources.
lambda = [400:10:700]*1e-9;
R = loadspectrum(lambda, 'redbrick');

sun = loadspectrum(lambda, 'solar');
lamp = blackbody(lambda, 2600);

xy_sun = lambda2xy(lambda, sun .* R)
xy_lamp = lambda2xy(lambda, lamp .* R)

%% 10.3.6 Gamma - making linear transition of brightness
wedge = [0:0.1:1];
figure;idisp(wedge)
figure;idisp( wedge .^ (1/2.2) ) % gamma can make the transition more linear


%% 10.4 Application Color Image

%% 10.4.1 Comparing Color Spaces

lambda = [400:5:700]*1e-9;
sun_ground = loadspectrum(lambda, 'solar');
figure;plot(lambda, sun_ground)

lambda = [400:5:700]*1e-9;
macbeth = loadspectrum(lambda, 'macbeth');
d65 = loadspectrum(lambda, 'D65') * 3e9;

clear Lab

% The greta macbeth color board has 24 colors, we can compare color spaces
% using different color methods.

 for i=1:18
 L = macbeth(:,i) .* d65;
 tristim = max(cmfrgb(lambda, L), 0);
 RGB = igamm(tristim, 0.45);

 XYZ(i,:) = colorspace('XYZ<-RGB', RGB);
 Lab(i,:) = colorspace('Lab<-RGB', RGB);
 end

xy = XYZ(:,1:2) ./ (sum(XYZ,2)*[1 1]);
ab = Lab(:,2:3);

figure;showcolorspace(xy', 'xy');
figure;showcolorspace(ab', 'ab'); % not sure why this doesn't work. i tried troubleshooting.
% corrections need to be done in the functions. I'll check version updates
% to avoid user erros.

%% 10.4.2 Shadow Removal
% So from my understanding shadows can have color is another light source
% is shining on them. For example, outside on a sunny day, a tree and grass
% is luminated by the sun and the sky but the shadow of the tree on the
% grass is only iluminatedby the sky. The sky has a different temperature
% (slightly blue) than the sun. Black bodies can be used to correct.
im1 = iread('parks.jpg');
idisp(im1)
im = iread('parks.jpg', 'gamma', 'sRGB');
figure;idisp(im)
%gs = invariant(im, 0.7, 'noexp');
% gs = RGB2IlluminationInvariant(im, 0.7, 'noexp');
gs = RGB2IlluminationInvariant(im1, 10000);

figure;idisp(gs)
%theta = esttheta(im)





