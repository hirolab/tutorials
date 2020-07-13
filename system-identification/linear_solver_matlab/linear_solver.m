%clc, clear all, close all
ts = 1/183
t = 0:ts:3;
x = sin(pi*2*3*t)-sin(pi*2*1*t)+2*sin(pi*2*0.5*t);
dx = diff(x)/ts
ddx = diff(dx)/ts
I = 1
B = 10
K = 100
C = 50
G = 3
F = I.*ddx+B.*dx(1:end-1)+K.*x(1:end-2)+C*sign(dx(1:end-1))+G

figure(1);plot(F/1000,'g');hold on; plot(x)
figure(2);plot(x(1:end-2),F);

A = [x(1:end-2)',dx(1:end-1)',ddx',sign(dx(1:end-1))',ones(length(ddx),1)]
Q = A\F'
