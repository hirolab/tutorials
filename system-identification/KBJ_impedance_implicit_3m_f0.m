classdef KBJ_impedance_implicit_3m_f0 < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        A, B
        Q, R, Ru, Ry
        
        stiff, damp, mass, mass2
        
        Deriv_blk
        nord = 5
    end
    
    methods
        
        function [obj, u, y, x0, P0, x] = KBJ_impedance_implicit_3m_f0(stiff, damp, mass, mass2, Fs, mocapCov, forceCov, ...
                                                                    xp, xf, xs, f, resCov, f0)

            n = obj.nord;
            
            % Build dataset ==================================
            dlen = size(xp, 1);
            
            rows = (10 : dlen-10)';
            u = [xp(rows + (n+1)/2), xf(rows + (n+1)/2), xs(rows + (n+1)/2)];
            uCov = diag([mocapCov, mocapCov, mocapCov]);
            y = [f(rows)];  % assume implicit function
            yCov = diag([forceCov, resCov]);

            x = zeros(length(rows), n*3);
            for i = 1 : n
                tshift = (n+1)/2 - i;
                x(:, i+[0, n, 2*n]) = [xp(rows+tshift), xf(rows+tshift), xs(rows+tshift)];
            end

            p0 = eye(n) * mocapCov;
            
            % Build system ===================================
                     
            a = [zeros(1, n); eye(n-1), zeros(n-1,1)];
            b = [1; zeros(n-1,1)];
                     
            time = ((n-1)/2 : -1 : -(n-1)/2)' / Fs;
            deriv_blk = inv((time.^(0:(n-1)) ./ factorial(0:(n-1))));
            deriv_blk = deriv_blk(1:3,:);
            
            x0 = x(1,:)';
            P0 = blkdiag(p0, p0, p0);         
            obj.A = blkdiag(a, a, a);
            obj.B = blkdiag(b, b, b);
            obj.Deriv_blk = blkdiag(deriv_blk, deriv_blk, deriv_blk(1:2,:));
            
            obj.Q = obj.B * uCov * obj.B';
            obj.R = yCov;  % force, imu
            obj.Ru = uCov;  % plate, foot, shank
            obj.Ry = yCov;
            
            obj.stiff = stiff;
            obj.damp = damp;
            obj.mass = mass;
            obj.mass2 = mass2;
            
            % Add ankle torque as state =====================
            x = [x, f0 * ones(length(rows),1), zeros(length(rows),1)];
            obj.A = blkdiag(obj.A, [1, 1; 0, 1]);
            obj.B = [obj.B; zeros(2, 3)];
%             obj.Deriv_blk = [obj.Deriv_blk, zeros(size(obj.Deriv_blk,1), 2)];
            obj.Q = blkdiag(obj.Q, zeros(2));
            
            x0 = x(1,:)';
            P0 = blkdiag(P0, eye(2));
        end
        
        function x2 = StateTransition(obj, x, u)
            if size(x, 2) == 1
                x2 = obj.A*x + obj.B*u;
            else
                x2 = x * obj.A' + u * obj.B';
            end
        end
        
        function h = Measurement(obj, x, v, uy)
            
            z = obj.Deriv_blk * x(1 : 3*obj.nord);
            xp = z(1:3);
            xf = z(4:6);
            xs = z(7:8);
            
            f0 = x(end-1);
            f = uy(end) + v(1);
            
            params = [obj.stiff, obj.damp, obj.mass, obj.mass2];
            res = KBJ_impedance_implicit_3m_f0.residual_fcn(...
                params, xp(1), xp(2), xp(3), xf(1), xf(2), xf(3), xs(1), xs(2), f, f0);

            h = res + v(2);
        end
        
        function dfdx = StateTransitionJacobian(obj, x, u)
            dfdx = obj.A;
        end
        
        function [dhdx, dhdv] = MeasurementJacobian(obj, x, v, uy)
            
            z = obj.Deriv_blk * x(1 : 3*obj.nord);
            xp = z(1:3);
            xf = z(4:6);
            xs = z(7:8);
            
            f0 = x(end-1);
            f = uy(end) + v(1);
            
            params = [obj.stiff, obj.damp, obj.mass, obj.mass2];
            [~, ~, res_z, res_f, res_f0] = KBJ_impedance_implicit_3m_f0.residual_fcn(...
                params, xp(1), xp(2), xp(3), xf(1), xf(2), xf(3), xs(1), xs(2), f, f0);
             
            dhdx = [res_z * obj.Deriv_blk, res_f0, 0];
            
            dhdv = [res_f, 1];
        end
        
        function [Sinv, Q, func] = wrap_cost_fcn(obj, X, P, U, Y)
            
            Q = [X(:,1:3*obj.nord) * obj.Deriv_blk', Y(:,1), X(:,end-1)];
            func = @(params, Q) KBJ_impedance_implicit_3m_f0.residual_fcn(params, ...
                Q(:,1), Q(:,2), Q(:,3), Q(:,4), Q(:,5), Q(:,6), Q(:,7), Q(:,8), Q(:,9), Q(:,10));
            
            nh = 1;
            ny = size(Y, 2);
            N = size(X, 1);
            Sinv = zeros(nh, nh, N);
            for i = 1 : N
                [dhdx, dhdv] = MeasurementJacobian(obj, X(i,:)', zeros(ny,1), [U(i,:), Y(i,:)]');
                S = dhdx * P(:,:,i) * dhdx' + dhdv * obj.Ry * dhdv';
                Sinv(:,:,i) = inv(S(1:nh,1:nh));
            end
        end
    end
    
    methods (Static)
        
        function [res, res_par, res_q, res_f, res_f0] = residual_fcn(params, sp, vp, ap, sf, vf, af, ss, vs, f, f0)
            % x: [sp, vp, ap, sf, vf, af, ss, vs]
            stiff = params(1);
            damp = params(2);
            mass = params(3);
            mass2 = params(4);
            
            res = mass*af + mass2*ap + damp*(vf - vs) + stiff*(sf - ss) - f + f0;
            res_par = [sf - ss, vf - vs, af, ap];
            res_q = [0, 0, mass2, stiff, damp, mass, -stiff, -damp];
            res_f = -1;
            res_f0 = 1;
        end
    end
    

end

