function [cost, grad] = cost_fcn(params, Rinv, X)

alf = params(1:2);
V1 = reshape(params(3:end), [], 1);
dlen = size(V1, 1);

x1 = X(:,1) - V1;
x2 = X(:,2);

V2 = -(alf(1)*x1.^2 + alf(2) - x2);
error = [V1, V2];

cost = 0;
gradV = zeros(dlen, 1);
gradAlf = [0, 0];
for i = 1 : dlen
    cost = cost + error(i,:) * Rinv * error(i,:)' / dlen;
    gradCost = 2 * error(i,:) * Rinv / dlen;
    gradAlf = gradAlf + gradCost * [0, 0; -x1(i)^2, -1];
    gradV(i) = gradCost * [1; 2*alf(1)*x1(i)];
end

grad = [gradAlf, reshape(gradV, 1, [])];

