classdef SimpleDynamicModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = untitled(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods (Static)
        function Xkk = state_fcn(X)
            % spring mass damper
            
%             A = [ 0.953488372093023, 0.009302325581395
%                  -9.302325581395349, 0.860465116279070];   % spring-mass-damper 100 1 0.1
             
            dt = 1/100;
            A = [   1,  dt, dt^2
                    0,  1,  dt
                    0,  0,  1];  % constant acceleration
            Xkk = A * X;
        end
        
        function Y = measurement_fcn(X)
            % read position
            
%             C = [1; 0];
%             C = [0.976744186046512, 0.004651162790698];
            C = [1, 0, 0; -1, 0, 0; 0.5, 0, 0; -0.3, 0, 0];
            Y = C * X;
        end
        

    end
end

