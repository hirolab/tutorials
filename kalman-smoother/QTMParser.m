classdef QTMParser
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Fs
        bodies
        markers
        
        pforce, ptorque
    end
    
    properties (SetAccess = private, Hidden = true)
        m  % measurement file from QTM

        R_ISB_Kistler = [1 0 0; 0 0 -1; 0 1 0]  % rotm to ISB from Kistler
    end
    
    methods
        function obj = QTMParser(pathMat)
            % Load mat file ============================
            d = load(pathMat);
            measName = fields(d);
            assert(length(measName) == 1, 'Multiple measurement files!')
            obj.m = d.(measName{1});
            obj.Fs = obj.m.FrameRate;
            
            for i = 1 : length(obj.m.RigidBodies.Name)
                bodies(i) = obj.loadBody(obj.m.RigidBodies.Name{i});
            end
            obj.bodies = bodies;
            
            for i = 1 : length(obj.m.Trajectories.Labeled.Labels)
                markers(i) = obj.loadMarker(obj.m.Trajectories.Labeled.Labels{i});
            end
            obj.markers = markers;
            
            
            FsP = obj.m.Force.Frequency;
            pforceP = obj.m.Force.Force' * obj.R_ISB_Kistler';
            ptorqueP = obj.m.Force.Moment' * obj.R_ISB_Kistler';
            
            inf_rows = any(~isfinite([pforceP, ptorqueP]), 2);
            pforceP(inf_rows,:) = interp1(find(~inf_rows), pforceP(~inf_rows,:), find(inf_rows), 'nearest', 'extrap');
            ptorqueP(inf_rows,:) = interp1(find(~inf_rows), ptorqueP(~inf_rows,:), find(inf_rows), 'nearest', 'extrap');

            ptimeP = (0:(size(pforceP,1)-1))' / FsP;
            time = (0:(obj.m.Frames-1))' / obj.Fs;

            [bf, af] = fir1(2000, (obj.Fs*0.4)/FsP);
            obj.pforce = interp1(ptimeP, filtfilt(bf, af, pforceP), time);
            obj.ptorque = interp1(ptimeP, filtfilt(bf, af, ptorqueP), time);
        end
    end
    
    
    methods (Access = private)
        function b = loadBody(obj, name)
            b.name = name;

            % extract Rigid Bodies
            [Lia, Locb] = ismember(b.name, obj.m.RigidBodies.Name);
            assert(all(Lia), 'Body not found')
            b.trans = permute(obj.m.RigidBodies.Positions(Locb,:,:) * 1e-3, [3 2 1]);
            b.quat = dcm2quat(reshape(obj.m.RigidBodies.Rotations(Locb,:,:), 3, 3, [])) .* [1, -1, -1, -1];
            % b.quat = rotm2quat(reshape(obj.m.RigidBodies.Rotations(Locb,:,:), 3, 3, []));
            b.residual = reshape(obj.m.RigidBodies.Residual(Locb,:,:), [], 1);

            % smooth quat
            quat = b.quat;
            for k = 2 : size(quat, 1)
                if sum(abs(quat(k,:) - quat(k-1,:))) > sum(abs(quat(k,:) + quat(k-1,:)))
                    quat(k:end,:) = -quat(k:end,:);
                end
            end
            b.quat = quat;
        end
        
        function b = loadMarker(obj, name)
            b.name = name;

            % extract Markers
            [Lia, Locb] = ismember(b.name, obj.m.Trajectories.Labeled.Labels);
            assert(all(Lia), 'Not all markers found')
            data = permute(obj.m.Trajectories.Labeled.Data(Locb,:,:), [3 2 1]);
            b.trans = data(:,1:3,:) * 1e-3;  % to [m]
            b.residual = data(:,4,:) * 1e-3;
        end
       
    end
    
    methods (Static)
        function plotBodyMarkers(b, m)
            nMarkers = numel(m);
            
            mtransAvg = zeros(nMarkers, 3);
            for k = 1 : nMarkers
                mtransAvg(k,:) = nanmean(quatrotate(b.quat, m(k).trans - b.trans));
            end
            tmp = randi(nMarkers, [10*nMarkers*nMarkers, 1]);  % lines inter-connecting all markers

            rotAvg = eye(3);
            transAvg = zeros(1, 3);

            plot3(mtransAvg(:,1), mtransAvg(:,2), mtransAvg(:,3), 'o')
            hold on, plot3(mtransAvg(tmp,1), mtransAvg(tmp,2), mtransAvg(tmp,3))
            plot3(transAvg(1), transAvg(2), transAvg(3), 'sk', 'markerfacecolor', 'k')
            text(transAvg(1), transAvg(2), transAvg(3), ['   ' b.name])
            quiver3(transAvg(1), transAvg(2), transAvg(3), rotAvg(1,1), rotAvg(1,2), rotAvg(1,3), 0.1, 'r')
            quiver3(transAvg(1), transAvg(2), transAvg(3), rotAvg(2,1), rotAvg(2,2), rotAvg(2,3), 0.1, 'g')
            quiver3(transAvg(1), transAvg(2), transAvg(3), rotAvg(3,1), rotAvg(3,2), rotAvg(3,3), 0.1, 'b')
            xlabel('x'), ylabel('y'), zlabel('z'), grid on, axis equal

            ds = 100;  % time-decimation
            for k = 1 : nMarkers
                mtrans = quatrotate(b.quat, m(k).trans - b.trans);
                mtrans0 = nanmean(mtrans);
                C = nancov(mtrans - mtrans0);
                if any(isnan(C(:)))  % marker or body never visible
                    warning('Marker not visible')
                    continue
                end
                [U,S,~] = svd(C);

                plot3(mtrans(1:ds:end,1), mtrans(1:ds:end,2), mtrans(1:ds:end,3), '.')
                hold on
                quiver3(mtrans0(1), mtrans0(2), mtrans0(3), U(1,1), U(2,1), U(3,1), 10*sqrt(S(1,1)), 'r')
                quiver3(mtrans0(1), mtrans0(2), mtrans0(3), U(1,2), U(2,2), U(3,2), 10*sqrt(S(2,2)), 'g')
                quiver3(mtrans0(1), mtrans0(2), mtrans0(3), U(1,3), U(2,3), U(3,3), 10*sqrt(S(3,3)), 'b')
                text(mtrans0(1), mtrans0(2), mtrans0(3), ['   ' m(k).name])
            end
        end
        
    end
end

