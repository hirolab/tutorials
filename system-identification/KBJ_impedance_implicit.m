classdef KBJ_impedance_implicit < handle
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        A, B, C, D, Jy
        Q, R
    end
    
    methods
        
        function obj = KBJ_impedance_implicit(A, B, C, D, Jy, Ry, Ru)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            obj.Jy = Jy;
            obj.Q = B * Ru * B';
            obj.R = blkdiag(Ru, Ry);
        end
        
        function x2 = StateTransition(obj, x, u)
            if size(x, 2) == 1
                x2 = obj.A*x + obj.B*u;
            else
                x2 = x * obj.A' + u * obj.B';
            end
        end
        
        function h = Measurement(obj, x, v, uy)
            if size(x, 2) == 1
                h = obj.C*x + [obj.D, obj.Jy] * (uy + v);
            else
                h = x * obj.C' + uy * [obj.D, obj.Jy]';
            end
        end
        
        function dfdx = StateTransitionJacobian(obj, x, u)
            dfdx = obj.A;
        end
        
        function [dhdx, dhdv] = MeasurementJacobian(obj, x, v, uy)
            dhdx = obj.C;
            dhdv = [obj.D, obj.Jy];
        end
    end
end

