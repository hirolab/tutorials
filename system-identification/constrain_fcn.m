function [c, ceq, gradC, gradCeq] = constrain_fcn(params, X)

dlen = size(X, 1);
alf = params(1:2);
V = reshape(params(3:end), [], 2);

x1 = X(:,1) - V(:,1);
x2 = X(:,2) - V(:,2);

c = [];
% % 0 == alf1 * (xv - v)**2 + alf2 - (yw - w)     ||    y = 1.0 * x.^2 + 1.0;
ceq = alf(1)*x1.^2 + alf(2) - x2;

if nargout > 2
    gradi = zeros(dlen, 2*dlen);
    for i = 1 : dlen
        gradi(i,i+[0, dlen]) = [-2*alf(1)*x1(i), 1];
    end

    gradC = [];
    gradCeq = [x1.^2, ones(dlen, 1), gradi]';
end