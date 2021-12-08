cd('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training')
clc; clear all; close all;
addpath('C:\Users\hirolab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\vision\mex')

%%
%cd('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training\TOOLBOX_calib')
addpath('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training\TOOLBOX_calib')
addpath('C:\Users\hirolab\Documents\MATLAB\Computer Vision Training\Calibration_Images')

%% IMPORTANT: If 'i'functions don't work run
% cd('C:\Users\hirolab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Machine Vision Toolbox for MATLAB\vision')
% make

%% 12.1 Obtaining an Image

%%  Working with unit8 
% when working with uint 8 arithmetic, data is clippled from 0 to 255
a = uint8(100)
b = uint8(200)
a+b
a-b

%% 12.1.1 image from files
street = iread('street.png'); %iread imports an image
about(street) %info of the image
imshow(street)     %pixes in matlab have an inverted notation to the book
street(200,300)
%%
streetd = iread('street.png', 'double'); %converted from uint(unsigned integer 8 byte, to float
about streetd

idisp(street) %idsp allows for interactive display of the image

%%

flowers = iread('flowers8.png');
about(flowers)
pix = flowers(276,318,:)
%idisp(flowers)

about pix
squeeze(pix)'

idisp( flowers(:,:,1) ) % red color plane index (white is high, black is absent)
%% Not Runnable but still useful functions
% seq = iread('seq/*.png');  %loads a sequence of wildcard name files 
% about(seq)
% loads nine images in PNG format from the folder seq. The result is a H×W×N matrix
% and the last index represents the image number within the sequence. That is
% seq(:,:,k) is the kth image in the sequence and is a 512 × 512 greyscale image

%% md = metadata
[im,md]=iread('street1.jpg');
md
%md.DigitalCamera

%% 12.1.2 Images from attached cameras

VideoCamera('?')
%%
% https://www.mathworks.com/help/supportpkg/usbwebcams/ug/acquire-images-from-webcams.html
webcamlist

cam = webcam(1)
%%
cam.AvailableResolutions
cam.Resolution = '640x480';
preview(cam)
%closePreview(cam)

%% These codes only work external cameras I will return to these
cam.size
im = cam.grab();
%% Movies
cam = Movie('traffic_sequence.mpg');
cam.size()
im = cam.grab();
about(im)
while 1
 im = cam.grab;
 if isempty(im) break; end
 image(im); drawnow
 end
%% Images from web

cam = AxisWebCamera('http://wc2.dartmouth.edu');
im = cam.grab();
idisp(im)

%% Images from map
ev = EarthView('key', YOUR_KEY); % need API key

%% Test patterns
im = testpattern('rampx', 256, 2);
idisp(im)
im = testpattern('siny', 256, 2);
idisp(im)
im = testpattern('squares', 256, 50, 25);
idisp(im)
im = testpattern('dots', 256, 256, 100);
idisp(im)

%%  Canvas
canvas = zeros(1000, 1000); % black pixel canvas
sq1 = 0.5 * ones(150, 150);
sq2 = 0.9 * ones(80, 80);

canvas = ipaste(canvas, sq1, [100 100]);
canvas = ipaste(canvas, sq2, [300 300]);
circle = 0.6 * kcircle(120);
canvas = ipaste(canvas, circle, [600, 200]);
canvas = iline(canvas, [100 100], [800 800], 0.8);

idisp(canvas)
%%
church = iread('church.png', 'grey');
figure;idisp(church)
figure;ihist( church )
% Histogram is centered, if too much exposure occured histogram would shift
% right, if under exposed histogram would shift left.
[n,v] = ihist(church);
[~,x] = peak(n, v)
about x
[~,x] = peak(n, v, 'scale', 25)

%% 12.3 monadic image-processing operations
% functions for each pixel
imd = idouble(church); % change uint8 to double precision 
im = iint(imd); % vice versa
flowers = iread('flowers8.png');
grey = imono(flowers); %grey scale
figure;idisp(grey)
color = icolor(grey); %return 3 dimension vector
%idisp(color)
color = icolor(grey, [1 0 0]); %just red colors
figure;idisp(color)

%thresholding
bright = (church >= 180);
figure;idisp(bright)

im = istretch(church);
im = inormhist(church); % forms of mapping and stretching histogram to use full range of pixels
figure;idisp(im)
figure;ihist(church, 'cdf')
im = igamm(church, 1/0.45);
figure;idisp(im)
im = igamm(church, 'sRGB');
figure;idisp(im)
idisp( church/64 )

%% Review of 10.3.6 Gamma before continuing with 12.4
wedge = [0:0.1:1];
figure;idisp(wedge)
figure;idisp( wedge .^ (1/2.2) ) % gamma can make the transition more linear

%% 12.4 Diadic Image Processing 
% where two input images result in a single output image
subject = iread('greenscreen.jpg', 'double');
idisp(subject)
linear = igamm(subject, 'sRGB');
idisp(linear)
[r,g] = tristim2cc(linear);

ihist(g)
mask = g < 0.45;
idisp(mask)
mask3 = icolor( idouble(mask) );
idisp(mask3 .* subject);
bg = iread('road.png', 'double');
bg = isamesize(bg, subject);
idisp(bg .* (1-mask3))

idisp( subject.*mask3 + bg.*(1-mask3) );
ipixswitch(mask, subject, bg);
%% Similar background removal but without green screen
vid = Movie('traffic_sequence.mpg', 'grey', 'double');
implay('traffic_sequence.mpg');

bg = vid.grab();

sigma = 0.02;
while 1
 im = vid.grab;
 if isempty(im) break; end; % end of fi le?
 d = im-bg;
 d = max( min(d, sigma), -sigma); % apply c(.)
 bg = bg + d;
 idisp(bg);
end

%% 12.5

%% Convolution 
%O = iconvolve(K, I); 
CC = 11
CC = 3
CC = 21
K = ones(CC,CC) / CC^2;
mona = iread('monalisa.png', 'double', 'grey');
idisp(mona)
idisp( iconvolve(mona, K) );

K = kgauss(5);
idisp( iconvolve(mona, K) );

about(K)
idisp( ismooth(mona, 5) )  % although blurring seems useless, it is actually an important process in image processing!
idisp( K );
surfl (-15:15, -15:15, K);
K = kcircle(8, 15);

%% edge detection
castle = iread('castle.png', 'double', 'grey');
idisp(castle)
p = castle(360,:);
plot(p)
plot(diff(p))  % differential in change of color is great indicator of edges (better than pixel color)

% convolution using the provided K is the opposite of smoothing.
K = [0.5 0 -0.5]; % colors can be changed with the sign of the slope K
idisp( iconvolve(castle, K), 'invsigned')

Du = ksobel
idisp( iconvolve(castle, Du), 'invsigned') % draws horizontal lines for edges
idisp( iconvolve(castle, Du'), 'invsigned') % draws vertical lines to find horizontal edges

Gu = iconvolve( Du, kgauss(sigma) , 'full'); %derivative of Gaussian
idisp(Gu)

Iu = iconvolve( castle, kdgauss(2) );idisp(Iu)
Iv = iconvolve( castle, kdgauss(2)' );idisp(Iv)
m = sqrt( Iu.^2 + Iv.^2 ); % magnitude of gradient at each pixel


figure; %for the sparsed quiver plot
th = atan2( Iv, Iu);
quiver(1:20:numcols(th), 1:20:numrows(th), ...
Iu(1:20:end,1:20:end), Iv(1:20:end,1:20:end))


%% smoothing the image before finding edges can facilitate process
castle = iread('castle.png', 'double', 'grey');
idisp( ismooth(castle, 5) )
Du = ksobel
idisp( iconvolve(ismooth(castle, 5), Du), 'invsigned') % draws horizontal lines for edges


%%

edges = icanny(castle, 2);
idisp(edges)
L = klaplace() % laplacian operator (is isotropic, 
% meaning it responds equally to any direction of edge

lap = iconvolve( castle, klog(2) );
idisp(lap,'signed')

p = lap(360,570:600);
plot(570:600, p, '-o');
zc = zcross(lap);

%% Image template (matching)

mona = iread('monalisa.png', 'double', 'grey');
T = mona(170:220, 245:295);
% perfect cases
sad(T, T) % sum absolute difference (SAD)
ssd(T, T) % sum square difference (SSD)
ncc(T, T) % normalized cross correlation

%one is darker
sad(T, T*0.9)
ssd(T, T* .9)

ncc(T, T*0.9) % NCC is unafected by brightness

% with an offset this time
zsad(T, T+0.1)
zssd(T, T+0.1)
zncc(T, T+0.1) %zero mean normalized cross correlation

% these methods do not work if the image changes in scale or rotation!!!!
%% where's waldo example 

crowd = iread('wheres-wally.png', 'double');
idisp(crowd)
wally = iread('wally.png', 'double');
idisp(wally); title('what he looks like')

S = isimilarity(wally, crowd, @zncc); % with this we find numerous candidates
idisp(S, 'colormap', 'jet', 'bar') 


% install: MATLAB Support for MinGW-w64 C/C++ Compiler 
% then go to the folder and run make (this will install iwindow)

[mx,p] = peak2(S, 1, 'npeaks', 5); % uses 2-D peak findin to compare candidates %
idisp(crowd);
plot_circle(p, 30, 'edgecolor', 'g')
plot_point(p, 'sequence', 'bold', 'textsize', 24, 'textcolor', 'y')

% number one is in fact waldo!
%% Nonlinear Functions
idisp(mona)
out = iwindow(mona, ones(7,7), 'var'); %computes variance of every pixel
mx = irank(mona, 1, 2);
idisp(mx)
med = irank(mona, 12, 2);
idisp(med)
% these procedures remove inpulse type noises. In this case, the cracks of
% mona lisa are gone

%% This method is better than blurring to elimiante impulse noise! Better image quality.
lena = iread('church.png','double','grey');
%lena = iread('lena.pgm', 'double'); couldnot find intended image
spotty = lena;
npix = prod(size(lena));
spotty(round(rand(5000,1)*(npix-1)+1)) = 0;
spotty(round(rand(5000,1)*(npix-1)+1)) = 1.0;
figure;idisp(spotty)
figure;idisp( irank(spotty, 5, 1) )

%% 12.6 Mathematical Morphology 
eg_morph1
idisp(im)
S = ones(5,5);
mn = imorph(im, S, 'min');
figure;morphdemo(im, S, 'min') %the new image is the minimum of the moving square window 5x5 
% only shapes that can contained the square survived

mx = imorph(mn, S, 'max');
idisp(mx) % errosion and dilation (changing size of square)

out = ierode(im, S);
out = idilate(im, S); % if any part is 1 then the whole section becomes 1
idisp(out)

%% 12.6.1 Noise Removal 

objects = iread('rice.png','double');
idisp(objects)

S = kcircle(3)
closed = iclose(objects, S);figure;idisp(closed)
clean = iopen(closed, S);figure;idisp(clean)
opened = iopen(objects, S);figure;idisp(opened) % best noise filter removal
closed = iclose(opened, S);figure;idisp(closed)


%% Boundary detection (using previous section 12.6.1)
eroded = imorph(clean, kcircle(1), 'min');
idisp(clean-eroded)

cleaneroded = clean-eroded

idisp(imcomplement(imbinarize(opened)))

image = imcomplement(imbinarize(opened))

out = hitormiss(cleaneroded, S);idisp(out)% didn't work

skeleton = ithin(image);
idisp(skeleton)

ends = iendpoint(skeleton);
joins = itriplepoint(skeleton);

%% 12.6.4 Distance Transform

im = testpattern('squares', 256, 256, 128);
idisp(im)

im = irotate(im, -3);idisp(im)
figure;
edges = icanny(im) > 0;
idisp(edges)

dx = distancexform(edges, 'euclidean');
dx = distancexform(edges);

%% 12.7 Shape Changing

mona = iread('monalisa.png');
[eyes,roi] = iroi(mona); %ROI  = region of interest
idisp(eyes) ; % interactive click and drag cursos to crop on command
roi

smile = iroi(mona, [265 342; 264 286]);
idisp(smile)

%% 12.7.2 Image resizing
roof = iread('roof.jpg', 'grey');
about(roof)

smaller = roof(1:7:end,1:7:end); % reducing by a factor of 7
about(smaller) % the problem with this is anti-aliasing
idisp(roof)
figure;idisp(smaller)

smaller = idecimate(roof, 7); %applies blurr before reducing image to avoid aliasing
figure;idisp(smaller)

% now bigger from the smaller file
bigger = ireplicate( smaller, 7 );
figure;idisp(bigger) % looks blocky
smoother = ismooth( bigger, 4);
idisp(smoother)% smoother

smaller = iscale(roof, 0.1);
bigger = iscale(smaller, 10);idisp(bigger)

%% Image Pyramid
p = ipyramid( imono(mona) )
idisp(p)

%% 12.7.4 Image Warping
% shrinking
mona = iread('monalisa.png', 'double', 'grey');
[Ui,Vi] = imeshgrid(mona);     %integer coordinates
[Up,Vp] = imeshgrid(400, 400); %integer coordinates
idisp(mona)

U = 4*(Up-100);  % float coordinates
V = 4*(Vp-200);
little_mona = interp2(Ui, Vi, mona, U, V);
figure;idisp(little_mona)

%rotating (using rotation matric math

mona = iread('monalisa.png', 'double', 'grey');
[Ui,Vi] = imeshgrid(mona);     %integer coordinates
[Up,Vp] = imeshgrid(600, 600); %integer coordinates
idisp(mona)

R = SO2(pi/6).R;
uc = 256; 
vc = 256;
U = R(1,1)*(Up-uc) + R(2,1)*(Vp-vc) + uc;
V = R(1,2)*(Up-uc) + R(2,2)*(Vp-vc) + vc;
twisted_mona = interp2(Ui, Vi, mona, U, V);
idisp(twisted_mona)
%% irotate function does this for us
twisted_mona = irotate(mona, rad2deg(pi/3),[]); % this one uses degrees 
idisp(twisted_mona)

%% Warping to undistort camera calibration 12.7.4

return
distorted = iread('Image18.tif', 'double');
[Ui,Vi] = imeshgrid(distorted);
Up = Ui; Vp = Vi;

%% http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example.html
% run:    calib_gui  first 
% write eimage names
% write image format
% click extract corners
% wintx = 5

load Bouguet

k = kc([1 2 5]); p = kc([3 4]);
u0 = cc(1); v0 = cc(2);
fpix_u = fc(1); fpix_v = fc(2);

u = (Up-u0) / fpix_u;
v = (Vp-v0) / fpix_v;


