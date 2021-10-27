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


% This Chapter was not very related to gait eyes, so I did not delve into the advance topics. 
% However, it is great for exploring relationships between color and
% energy or modifying colors. Gait eyes doesn't use colors (for now). So
% I'll stop here. 
