classdef KBJ_impedance_implicit_3m_aug < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        A, B
        Q, R, Ru, Ry
        T
        
        stiff, damp, mass, mass2
        
        Deriv_blk
    end
    
    methods
        
        function [obj, u, y, x0, P0, x] = KBJ_impedance_implicit_3m_aug(stiff, damp, mass, mass2, Fs, mocapCov, forceCov, imuCov, ...
                                                                    xp, xf, xs, f, xpdd)

            n = 5;
            nu = 3;
            T = eye(n) / 1e-2;
            obj.T = blkdiag(T, T, T, diag(1./[stiff, damp, mass, mass2]));
%             obj.T = eye(n*3 + 4);
            
            % Build dataset ==================================
            dlen = size(xp, 1);
            
            rows = (10 : dlen-10)';
            u = [xp(rows + (n+1)/2), xf(rows + (n+1)/2), xs(rows + (n+1)/2)];
            uCov = diag([mocapCov, mocapCov, mocapCov]);
            y = [f(rows), xpdd(rows)];  % assume implicit function
            yCov = diag([forceCov, imuCov]);

            x = zeros(length(rows), n*3);
            for i = 1 : n
                tshift = (n+1)/2 - i;
                x(:, i+[0, n, 2*n]) = [xp(rows+tshift), xf(rows+tshift), xs(rows+tshift)];
            end
            % add the params

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
%             obj.R(1) = obj.R(1) + 1e-2;
            obj.Ru = uCov;  % plate, foot, shank
            obj.Ry = yCov;
            
            obj.stiff = stiff;
            obj.damp = damp;
            obj.mass = mass;
            obj.mass2 = mass2;
            
            % augment with impedance params
            x = [x, repmat([stiff, damp, mass, mass2], [length(rows), 1])];
            x0 = x(1,:)';
            P0 = blkdiag(P0, diag([stiff*1e-4, damp*1e-4, 1*mass, 1e-4*mass2].^2));
            obj.A = blkdiag(obj.A, eye(4));
            obj.B = [obj.B; zeros(4, nu)];
            obj.Q = blkdiag(obj.Q, 0*diag([stiff, damp, mass, mass2].^2));
            
            % transform system state
            x = x * obj.T';
            x0 = obj.T * x0;
            P0 = obj.T * P0 * obj.T';
            obj.A = obj.T * obj.A / obj.T;
            obj.B = obj.T * obj.B;
            obj.Q = obj.T * obj.Q * obj.T';
        end
        
        function x2 = StateTransition(obj, x, u)
            if size(x, 2) == 1
                x2 = obj.A*x + obj.B*u;
            else
                x2 = x * obj.A' + u * obj.B';
            end
        end
        
        function h = Measurement(obj, x, v, uy)
            
            x = obj.T \ x;
            z = obj.Deriv_blk * x(1:end-4);
%             z = obj.Deriv_blk * x;
            xp = z(1:3);
            xf = z(4:6);
            xs = z(7:8);
            
            f = uy(end-1) + v(1);
            xpdd = uy(end) + v(2);
            
            params = x(end-3:end)';
%             params = [obj.stiff, obj.damp, obj.mass, obj.mass2];
            res = KBJ_impedance_implicit_3m_aug.residual_fcn(...
                params, xp(1), xp(2), xp(3), xf(1), xf(2), xf(3), xs(1), xs(2), f);

            h = [res
                 xp(3) - xpdd];
        end
        
        function dfdx = StateTransitionJacobian(obj, x, u)
            dfdx = obj.A;
        end
        
        function [dhdx, dhdv] = MeasurementJacobian(obj, x, v, uy)
            
            x = obj.T \ x;
            nh = 2;
            z = obj.Deriv_blk * x(1:end-4);
%             z = obj.Deriv_blk * x;
            xp = z(1:3);
            xf = z(4:6);
            xs = z(7:8);
            
            f = uy(end-1) + v(1);
            
            params = x(end-3:end)';
%             params = [obj.stiff, obj.damp, obj.mass, obj.mass2];
            [~, res_par, res_z, res_f] = KBJ_impedance_implicit_3m_aug.residual_fcn(...
                params, xp(1), xp(2), xp(3), xf(1), xf(2), xf(3), xs(1), xs(2), f);
             
            dhdx = [res_z
                 0, 0, 1, 0, 0, 0, 0, 0] * obj.Deriv_blk;
            dhdx = [dhdx, [res_par; zeros(1, 4)]];
            
            dhdv = [res_f, 0
                    0, -1];
                
            dhdx = dhdx / obj.T;
        end
        
        function [Sinv, Q, func] = wrap_cost_fcn(obj, X, P, U, Y)
            
            Q = [X * obj.Deriv_blk', Y(:,1)];
            func = @(params, Q) KBJ_impedance_implicit_3m_aug.residual_fcn(params, ...
                Q(:,1), Q(:,2), Q(:,3), Q(:,4), Q(:,5), Q(:,6), Q(:,7), Q(:,8), Q(:,9));
            
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
        
        function [res, res_par, res_q, res_f] = residual_fcn(params, sp, vp, ap, sf, vf, af, ss, vs, f)
            % x: [sp, vp, ap, sf, vf, af, ss, vs]
            stiff = params(1);
            damp = params(2);
            mass = params(3);
            mass2 = params(4);
            
            res = mass*af + mass2*ap + damp*(vf - vs) + stiff*(sf - ss) - f;
            res_par = [sf - ss, vf - vs, af, ap];
            res_q = [0, 0, mass2, stiff, damp, mass, -stiff, -damp];
            res_f = -1;
        end
    end
    

end

