function [Xs, Ps, X, Ppost, E] = smooth_LTI_implicit(A, B, C, D, Jy, Q, R, u, y, xk, P)

nx = size(A, 1);
ny = size(Jy, 1);
N = size(u, 1);
X = zeros(N, length(xk));  % [q, w, wd]
X(1,:) = xk;

Rsqrt = sqrtm(R);
Psqrt = sqrtm(P);
Qsqrt = sqrtm(Q);

E = zeros(N, ny);
Pprior = zeros(nx, nx, N);
Ppost = zeros(nx, nx, N);
Ypost = zeros(N, ny);
Ypre = zeros(N, ny);

Pprior_sqrt = zeros(nx, nx, N);
Ppost_sqrt = zeros(nx, nx, N);


for k = 1 : N  
    Pprior(:,:,k) = P;
    Pprior_sqrt(:,:,k) = Psqrt;
    
    % correct
    Ypre(k,:) = C * xk + D * u(k,:)' + Jy * y(k,:)';
%     Ypre(k,:) = C * xk + D * u(k,:)';
    e = - Ypre(k,:)';  % innovation via measurement data (will be weight-added to X)
%     e = y(k,:)' - Ypre(k,:)';  % innovation via measurement data (will be weight-added to X)
    S = C * P * C' + R;  % cov of e
       
    K = P * C' / S;
    xk = xk + K * e;  % t|t
%     P = P - K * S * K';  % lower norm(P - P') after many iters    
    [~, RR] = qr([Rsqrt, C*Psqrt; zeros(nx,ny), Psqrt]');
    Psqrt = RR(end-nx+1:end, end-nx+1:end)';
    P = real(Psqrt * Psqrt');
    
    X(k,:) = xk;
    E(k,:) = e;
    Ypost(k,:) = C * xk + D * u(k,:)' + Jy * y(k,:)';
%     Ypost(k,:) = C * xk + D * u(k,:)';
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
    xk = X(k,:)' + J*(xk - (A*X(k,:)' + B*u(k,:)'));
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

