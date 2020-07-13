function cost = MLE_cost(params, hfunc, Sinv)
% Maximum Likelihood Estimation of impedance params
%   Given smoothed states and covariance of error, calculate 
%   Input:
%       hfunc: nonlinear residual function 
%       S: covariance of the residual

res = hfunc(params);
dlen = size(Sinv, 3);

cost = 0;
for i = 1 : dlen
    cost = cost + res(i,:) * Sinv(:,:,i) * res(i,:)' / dlen;
end


