function [cost, jac] = MLE_cost(params, hfunc, Sinv, Q)
% Maximum Likelihood Estimation of impedance params
%   Given smoothed states and covariance of error, calculate 
%   Input:
%       hfunc: nonlinear residual function 
%       S: covariance of the residual

dlen = size(Sinv, 3);

cost = 0;
jac = zeros(1, length(params));
for i = 1 : dlen
    [res, grad] = hfunc(params, Q(i,:));
    cost = cost + res' * Sinv(:,:,i) * res / dlen;
    jac = jac + 2 * res' * Sinv(:,:,i) * grad / dlen;
end

