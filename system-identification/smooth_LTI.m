function [Xs, Ps, X, Ppost, E] = smooth_LTI(A, B, C, D, Q, R, u, y, xk, P)

nx = size(A, 1);
N = size(u, 1);
X = zeros(N, length(xk));  % [q, w, wd]
X(1,:) = xk;

E = zeros(N, 1);
Pprior = zeros(nx, nx, N);
Ppost = zeros(nx, nx, N);
Ypost = zeros(N, 1);
Ypre = zeros(N, 1);

for k = 1 : N  
    Pprior(:,:,k) = P;
    
    % correct
    Ypre(k,:) = C * xk + D * u(k,:)';
    e = y(k,:)' - Ypre(k,:)';  % innovation via measurement data (will be weight-added to X)
    S = C * P * C' + R;  % cov of e
    
    K = P * C' / S;
    xk = xk + K * e;  % t|t
    P = P - K * S * K';  % lower norm(P - P') after many iters
    
    X(k,:) = xk;
    E(k,:) = e;
    Ypost(k,:) = C * xk + D * u(k,:)';
    Ppost(:,:,k) = P;
    
    % update
    xk = A * xk + B * u(k,:)';  % t+1|t
    P = A * P * A' + Q;
end

% Pprior(:,:,N+1) = P;
P = Ppost(:,:,N);
xk = X(N,:)';

M = zeros(nx, nx, N);
M(:,:,N) = (eye(nx) - K*C)*A*Ppost(:,:,N-1);

Xs = zeros(N, nx);
Ps = zeros(nx, nx, N);
Xs(end,:) = X(end,:);
Ps(:,:,end) = Ppost(:,:,end);
% backward recursion
for k = N-1 : -1 : 1
    J = Ppost(:,:,k) * A' / Pprior(:,:,k+1);
    xk = X(k,:)' + J*(xk - A*X(k,:)' - B*u(k,:)');
    P = Ppost(:,:,k) + J*(P - Pprior(:,:,k+1)) * J';
    
    if k > 1
        Jm = Ppost(:,:,k-1) * A' / Pprior(:,:,k); % t-1
        M(:,:,k) = Ppost(:,:,k-1)*Jm' + J*(M(:,:,k+1) - A*Ppost(:,:,k))*Jm';
    end
    
    Xs(k,:) = xk;
    Ps(:,:,k) = P;
end


