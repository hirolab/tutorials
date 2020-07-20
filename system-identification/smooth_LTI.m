function [Xs, Ps, X, Ppost] = smooth_LTI(A, B, C, D, Q, R, u, y, xk, P)

nx = size(A, 1);
ny = size(y, 2);
N = size(u, 1);
X = zeros(N, length(xk));  % [q, w, wd]
X(1,:) = xk;

Rsqrt = sqrtm(R);
Psqrt = sqrtm(P);
Qsqrt = sqrtm(Q);

Pprior = zeros(nx, nx, N);
Ppost = zeros(nx, nx, N);

Pprior_sqrt = zeros(nx, nx, N);
Ppost_sqrt = zeros(nx, nx, N);

for k = 1 : N-1 
    
    % update
    xk = A * xk + B * u(k,:)';  % t+1|t
    P = A * P * A' + Q;
    
    Pprior(:,:,k) = P;
    Pprior_sqrt(:,:,k) = Psqrt;
    
    % correct
    ypre = C * xk + D * u(k,:)';
    e = y(k,:)' - ypre;  % innovation via measurement data (will be weight-added to X)
    S = C * P * C' + R;  % cov of e
       
    K = P * C' / S;
    xk = xk + K * e;  % t|t
    P = P - K * C * P;  % lower norm(P - P') after many iters    
%     [~, RR] = qr([Rsqrt, C*Psqrt; zeros(nx,ny), Psqrt]');
%     Psqrt = RR(end-nx+1:end, end-nx+1:end)';
%     P = real(Psqrt * Psqrt');
    
    X(k,:) = xk;
    Ppost(:,:,k) = P;
    Ppost_sqrt(:,:,k) = Psqrt;
    

%     [~, RR] = qr([Psqrt'*A'; Qsqrt']);
%     Psqrt = RR(1:nx,:)';
%     P = real(Psqrt * Psqrt');
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

return

Xs(end,:) = X(end,:);
Ps(:,:,end) = Ppost(:,:,end);
% backward recursion
for k = N-1 : -1 : 1
    J = Ppost(:,:,k) * A' / Pprior(:,:,k+1);
    xk = X(k,:)' + J*(xk - (A*X(k,:)' + B*u(k,:)'));
%     P = Ppost(:,:,k) + J*(P - Pprior(:,:,k+1)) * J';
    [~, RR] = qr([Ppost_sqrt(:,:,k)' * A', Ppost_sqrt(:,:,k)'
                  Qsqrt',                  zeros(nx,nx)
                  zeros(nx,nx),            Psqrt'*J']);
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

