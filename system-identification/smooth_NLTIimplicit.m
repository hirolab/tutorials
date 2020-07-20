function [Xs, Ps, Xf, Pf] = smooth_NLTIimplicit(obj, u, y, x0, P0)
% function [Xs, Ps, X, Ppost, E] = smooth_NLTIimplicit(A, B, Q, hfunc, Ru, Ry, u, y, xk, P)
% hfunc: residual function [hk, dh/dx, dh/du, dh/dy, dh/dparam] <- @(x, u, y)

% A = obj.A;  % assume invariant
Q = obj.Q;
N = size(u, 1);
nx = length(x0);
ny = size(y, 2);

% ukf = unscentedKalmanFilter(...
%     @obj.StateTransition,... % State transition function
%     @obj.Measurement,... % Measurement function
%     x0,...
%     'MeasurementNoise', obj.R, ...
%     'ProcessNoise', obj.Q, ...
%     'StateCovariance', P0, ...
%     'HasAdditiveMeasurementNoise', false, ...
%     'HasAdditiveProcessNoise', true);

ukf = extendedKalmanFilter(...
    @obj.StateTransition,... % State transition function
    @obj.Measurement,... % Measurement function
    x0,...
    'MeasurementNoise', obj.R, ...
    'ProcessNoise', obj.Q, ...
    'StateCovariance', P0, ...
    'HasAdditiveMeasurementNoise', false, ...
    'HasAdditiveProcessNoise', true, ...
    'StateTransitionJacobianFcn', @obj.StateTransitionJacobian, ...
    'MeasurementJacobianFcn', @obj.MeasurementJacobian);
%{
    , ...
    'StateTransitionJacobianFcn', @obj.StateTransitionJacobian, ...
    'MeasurementJacobianFcn', @obj.MeasurementJacobian
    %}

Xf = zeros(N, nx);
Pf = zeros(nx, nx, N);
Xp = zeros(N, nx);
Pp = zeros(nx, nx, N);
Xs = zeros(N, nx);
Ps = zeros(nx, nx, N);

for k = 1 : N
    [Xf(k,:), Pf(:,:,k)] = correct(ukf, zeros(ny, 1), [u(k,:)'; y(k,:)']);
    [Xp(k,:), Pp(:,:,k)] = predict(ukf, u(k,:)');
end

xk = Xf(N,:)';
P = Pf(:,:,N);
Psqrt = real(sqrtm(P));
Qsqrt = real(Q);

Xs(end,:) = Xf(end,:);
Ps(:,:,end) = Pf(:,:,end);
% backward recursion
for k = N-1 : -1 : 1
    A = StateTransitionJacobian(obj, xk, u(k,:)');
    
    Ppost_sqrt = real(sqrtm(Pf(:,:,k)));
    J = Pf(:,:,k) * A' / Pp(:,:,k);
    xk = Xf(k,:)' + J*(xk - Xp(k,:)');
    [~, RR] = qr([Ppost_sqrt' * A', Ppost_sqrt'
                  Qsqrt',                  zeros(nx,nx)
                  zeros(nx,nx),            Psqrt'*J']);
    Psqrt = RR(nx + (1:nx), (nx+1):end)';
    P = real(Psqrt * Psqrt');
    
    Xs(k,:) = xk;
    Ps(:,:,k) = P;
end

return



nx = size(A, 1);
hk = hfunc(xk, u(1,:)', y(1,:)');
nh = length(hk);
N = size(u, 1);
X = zeros(N, length(xk));  % [q, w, wd]
X(1,:) = xk;

Psqrt = sqrtm(P);
Qsqrt = sqrtm(Q);

E = zeros(N, nh);
Pprior = zeros(nx, nx, N);
Ppost = zeros(nx, nx, N);

Pprior_sqrt = zeros(nx, nx, N);
Ppost_sqrt = zeros(nx, nx, N);


for k = 1 : N  
    Pprior(:,:,k) = P;
    Pprior_sqrt(:,:,k) = Psqrt;
    
    % correct
    [hk, C, D, J] = hfunc(xk, u(k,:)', y(k,:)');
    R = D * Ru * D' + J * Ry * J';
    e = - hk;  % innovation via measurement data (will be weight-added to X)
    S = C * P * C' + R;  % cov of e
       
    K = P * C' / S;
    xk = xk + K * e;  % t|t
%     P = P - K * S * K';  % lower norm(P - P') after many iters    
    [~, RR] = qr([sqrtm(R), C*Psqrt; zeros(nx,nh), Psqrt]');
    Psqrt = RR(end-nx+1:end, end-nx+1:end)';
    P = real(Psqrt * Psqrt');
    
    X(k,:) = xk;
    E(k,:) = e;
    Ppost(:,:,k) = P;
    Ppost_sqrt(:,:,k) = Psqrt;
    
    % update
    xk = A * xk + B * u(k,:)';  % t+1|t
%     P = A * P * A' + Q;
    [~, RR] = qr([Psqrt'*A'; Qsqrt']);
    Psqrt = RR(1:nx,:)';
    P = real(Psqrt * Psqrt');
end

% Pprior(:,:,N+1) = P;
P = Ppost(:,:,N);
Psqrt = real(sqrtm(P));
xk = X(N,:)';

% M = zeros(nx, nx, N);
% M(:,:,N) = (eye(nx) - K*C)*A*Ppost(:,:,N-1);

Xs = zeros(N, nx);
Ps = zeros(nx, nx, N);
Ps_sqrt = zeros(nx, nx, N);

Xs(end,:) = X(end,:);
Ps(:,:,end) = Ppost(:,:,end);
% backward recursion
for k = N-1 : -1 : 1
    J = Ppost(:,:,k) * A' / Pprior(:,:,k+1);
    xk = X(k,:)' + J*(xk - A*X(k,:)' - B*u(k,:)');
%     P = Ppost(:,:,k) + J*(P - Pprior(:,:,k+1)) * J';
    [~, RR] = qr([Ppost_sqrt(:,:,k)' * A', Ppost_sqrt(:,:,k)'; Qsqrt', zeros(nx,nx); zeros(nx,nx), Psqrt'*J']);
    Psqrt = RR(nx + (1:nx), nx+1:end)';
    P = real(Psqrt * Psqrt');
    
%     if k > 1
%         Jm = Ppost(:,:,k-1) * A' / Pprior(:,:,k); % t-1
%         M(:,:,k) = Ppost(:,:,k-1)*Jm' + J*(M(:,:,k+1) - A*Ppost(:,:,k))*Jm';
%     end
    
    Xs(k,:) = xk;
    Ps(:,:,k) = P;
    Ps_sqrt(:,:,k) = Psqrt;
end


