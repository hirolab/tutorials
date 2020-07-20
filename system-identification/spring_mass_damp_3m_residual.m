function [hk, C, D, J, grad] = spring_mass_damp_3m_residual(xk, uk, yk, params)
% Spring-mass-damper in impedance form, with force and acceleration 
%   measurement in implicit form
%
% xk: buffer of [xp, xf, xs]
% uk: [xp, xf, xs] at k+(n+1)/2
% yk: [f, xpdd]
% params: [stiff, damp, mass]

stiff = params(1);
damp = params(2);
mass = params(3);
mass2 = params(4);

Fs = 183;
deriv_blk = [[     0,   0,    1,    0,     0]
             [ -1/12, 2/3,    0, -2/3,  1/12] * Fs
             [ -1/12, 4/3, -5/2,  4/3, -1/12] * Fs^2];
ord = size(deriv_blk, 2);

xp = deriv_blk * xk(1:ord);
xf = deriv_blk * xk((1:ord) + ord);
xs = deriv_blk(1:2,:) * xk((1:ord) + 2*ord);
f = yk(1);
xpdd = yk(2);

hk = [mass*xf(3) + mass2*xp(3) + damp*(xf(2) - xs(2)) + stiff*(xf(1) - xs(1)) - f
      xp(3) - xpdd];

dhdx = [0, 0, mass2, stiff, damp, mass, -stiff, -damp
        0, 0, 1, 0, 0, 0, 0, 0];
C = dhdx * blkdiag(deriv_blk, deriv_blk, deriv_blk(1:2,:));  % dh/dx
D = zeros(2, 3);  % dh/du
J = [1, 0; 0, 1];  % dh/dv, where v: [fv, app]

grad = [xf(1)-xs(1), xf(2)-xs(2), xf(3), xp(3)
        0, 0, 0, 0];  % dh/dparam
