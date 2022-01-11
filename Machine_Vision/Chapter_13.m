cd('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training')
clc; clear all; close all;
addpath('C:\Users\hirolab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\vision\mex')


%% make i-codes work
cd('C:\Users\hirolab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\vision')

addpath('C:\Users\hirolab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\vision')

make

%%VLfeat 
run('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training\src\vlfeat-devel\toolbox\vl_setup.m')

%% Image Feature Extraction

%% 13.1 Region Features

%% 13.111 Classification

%% Grey level classification
castle = iread('castle.png', 'double');
figure; idisp(castle);
figure;
idisp(castle >= 0.7)

ithresh(castle) % this threshold toolbox has an interactive slider

ithresh(castle*0.8) % different intensity day
ihist(castle);

%% Otsu method to determine optimum threshold
castle = iread('castle.png', 'double');
t = otsu(castle)
% however, this isnt alwasy perfect
castle = iread('castle2.png', 'double');
figure; idisp(castle);
t = otsu(castle)
idisp(castle >= 0.7)
ihist(castle);

%% one alternative to local highlights is to use local thresholds (Niblack algorithm)

 t = niblack(castle, -0.1, 30);
 figure;idisp(t)
 figure;idisp(castle >= t)
 
 %% Image Segmentation - install vl_mser to use following code!
 % run src/vlfeat-0.9.21/toolbox/vl_setup.m

 % to continue, download and install Vl-feat from: https://www.vlfeat.org/download.html
%  figure;
%  [mser,nsets] = imser(castle, 'area', [100 200]);
% idisp(mser, 'colormap', 'jet')

% will loop until all cases are completed 
% I added a counter to not lose patience
im = iread('castle_sign2.png', 'grey', 'double');
      [label,n] = imser(im, 'light');
      idisp(label)

%% 13.1.1.2 Color Classification

im_targets = iread('yellowtargets.png');
figure;idisp(im_targets)

im_garden = iread('tomato_124.jpg');
figure;idisp(im_garden)

%function colorkmeans fi rst maps each color pixel to a point on the
%xy- or a*b*-chromaticity plane. Then the k-means algorithm is used to fi nd clusters of
%points o

[cls, cab,resid] = colorkmeans(im_targets, 3,'spread'); % this line gave me trouble, I was tired of fixing codes by then.

%% 13.1.2 Representation

im = iread('multiblobs.png');
idisp(im)
[label, m] = ilabel(im);
idisp(label,'random')% labeling blobs
m  % number of blobs
idisp(label, 'colormap', jet, 'bar')
%% display only blob #3
reg3 = (label==3);
idisp(reg3) 

sum(reg3(:)) % number of pixels in the blob

%% connectivty analysi
[label, m, parents, cls] = ilabel(im);

parents
idisp(label, 'colormap', jet, 'bar')
targets_label = ilabel(targets_binary);
idisp(targets_label, 'colormap', 'jet');




%% 13.1.2.1 Graph based segmentation
cd('C:\Users\hirolab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\contrib\graphseg')
cd('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training\src\segment')
 cd('C:\Users\hirolab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\vision\mex')

im = iread('58060.jpg');
figure;idisp(im)
[label, m] = igraphseg(im, 1500, 100, 0.5);

%% Bounding boxes and moments
sharks = iread('sharks.png');
idisp(sharks)
[label, m] = ilabel(sharks);
blob = (label == 2);
figure;idisp(blob);
sum(blob(:))
[v,u] = find(blob);
about(u)

%bounds
umin = min(u)
umax = max(u)
vmin = min(v)
vmax = max(v)
plot_box(umin, vmin, umax, vmax, 'g')% look at figure

m00 = mpq(blob, 0, 0)
uc = mpq(blob, 1, 0) / m00
vc = mpq(blob, 0, 1) / m00
hold on; plot(uc, vc, 'gx', uc, vc, 'go');
u20 = upq(blob, 2, 0); u02 = upq(blob, 0, 2); u11 = upq(blob, 1, 1);
J = [ u20 u11; u11 u02]
plot_ellipse(4*J/m00, [uc, vc], 'b');
lambda = eig(J)

%% Blob features
%targets_binary = sharks;
f = imoments(blob)
f.uc
f.theta
f.aspect
f.moments.m00
f.moments.u11
fv = iblobs(sharks, 'boundary', 'class', 1);
idisp(sharks); hold on;
fv(1). plot_box('g')

fv(1)
about(fv(1).edge)

%% edgels (perimeter)
fv(1).edge(:,1:5)
idisp(sharks)
fv.plot_boundary('r')
fv.plot_centroid()
% Ciruclarity can also be used to identify objects!
[r,th] = fv(2).boundary; % perimeter and angle with respect to the centroid
plot([r th])
hold on
for f=fv
[r,t] = f.boundary();
plot(r/sum(r));
end
b = fv.boundary
RegionFeature.boundmatch(b(:,2), b) % similar blobs 
% blobs 2 is similae to 3 and 4

%% Character recognition
castle = iread('castle.png');
words = ocr(castle, [420 300 580 420])
words.Text

idisp(castle)
plot_box('matlab', words.WordBoundingBoxes, 'y')

%% Line Features
im = iread('5points.png', 'double');
idisp(im)

im = testpattern('squares', 256, 256, 128);
figure;idisp(im)

im = irotate(im, -10);
figure;idisp(im)
hold on;
edges = icanny(im);
figure;plot(edges)

h = Hough(edges)
lines = h.lines()
h = Hough(edges, 'suppress', 5)
lines = h.lines()
idisp(im);
h.plot('b')
%% lines (2)
im = iread('church.png', 'grey', 'double');
idisp(im);hold on;
edges = icanny(im);
h = Hough(edges, 'suppress', 10);
lines = h.lines();
h.plot('b')

idisp(im, 'dark');
lines(1:10).plot('g');

lines = lines.seglength(edges);
lines(1)
k = find( lines.length > 80);
lines(k).plot('b--')

%% 13.3  Points Features
b1 = iread('building2-1.png', 'grey', 'double');
idisp(b1)

C = icorner(b1, 'nfeat', 200)
figure
idisp(b1, 'dark'); %overlay corner points
C.plot('ws');
% to distribute points evenly and avoid clustering
figure;
Cs = icorner(b1, 'nfeat', 200, 'suppress', 10);
idisp(b1, 'dark'); %overlay corner points
Cs.plot('ws');

length(C)
C(1:4)

figure;
[C,strength] = icorner(b1, 'nfeat', 200);
idisp(strength, 'invsigned')

%% another point of view
b2 = iread('building2-2.png', 'grey', 'double');

C2 = icorner(b2, 'nfeat', 200);
idisp(b2,'dark')
C2.plot('ws');


%% Scale-space feature selection.
addpath("src")
im = iread('lena.png', 'double', 'grey');
idisp(im)

[G,L,s] = iscalespace(im, 50, 2);
idisp(L(:,:,5), 'invsigned')

%f = iscalemax(L, s)

[G,L] = iscalespace(im, 8, 8);
figure;idisp(G, 'flatten', 'wide', 'square');
figure;idisp(L, 'flatten', 'wide', 'square', 'invsigned');

%% Scale point feature
sf1 = isurf(b1, 'nfeat', 200)
sf1(1)
idisp(b1, 'dark');
sf1.plot_scale('g', 'clock')

figure;hist(sf1.scale, 100);










