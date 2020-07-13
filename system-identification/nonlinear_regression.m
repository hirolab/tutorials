clear, 

N = 400;
Fs = 183;
[bf, af] = fir1(10, 5/(Fs/2));
x1 = 2*filtfilt(bf, af, randn(N,1));
y = 1.0 * x1.^2 + 1.0;

% 0 == alf * (xv - v) - (yw - w)
Xn = [x1, y];
Rv = diag([2e-1, 2e-1].^2);
Rinv = inv(Rv);
V = mvnrnd([0; 0], Rv, N);
Xv = Xn + V;

sol1 = [Xv(:,1).^2, ones(N, 1)] \ Xv(:,2);
%
xn = [1; 1; V(:,1)];
x0 = [1; 1; zeros(N,1)];
lb = [0; 0; -ones(N,1)];
ub = [2; 2; ones(N,1)];
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'UseParallel', false, 'MaxFunctionEvaluations', 3e4, ...
    'SpecifyObjectiveGradient', true, ...
    'MaxIterations', 1e5, 'CheckGradients',false);

fun = @(x) cost_fcn(x, Rinv, Xv);
nonlcon = @(x) constrain_fcn(x, Xv);

tic
[x,fval,exitflag,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
toc

alf1 = x(1:2);
V1 = reshape(x(3:end), [], 1);

[sol1, alf1]
% [std(V1), std(V(:,1))]
% figure, plot([Xn(:,1), Xv(:,1), Xv(:,1)-V1])

%%


fun(xn)
[~, ceq] = nonlcon(xn)

