%% Model Ankle Dynamics
%% Model analytically
clear, clc

q = sym('q', [3, 1]);  % ankle angle
a = sym('a', [3, 1]);  % position of spring-damper on shank, shank frame
b = sym('b', [3, 1]);  % position of spring-damper on foot, foot frame

R = Ry(q(2)) * Rz(q(3)) * Rx(q(1));  % rotation shank to foot (YZX Euler)
% R = quat2rotm_sym([1, 0.5*q.']);  % small angle approximation
d = R * b - a;  % spring displacement
H = jacobian(d, q);  % jacobian to convert velocity/force between spring and ankle spaces
HH = reshape(jacobian(H(:), q), size(H,1), size(H,2), []);

syms stiff eabs  % spring and natural length of spring
e = eabs * d / sqrt(d.' * d);  % 
T = H.' * (stiff*(d - e));  % torque on ankle due to spring
K = jacobian(T, q);  % ankle stiffness due to spring
simplify(K - K.')  % equals zero for any parameter value

syms damp  % damping of spring
% B = -damp*(H.'*H);
qd = sym('qd', [3, 1]);  % ankle velocity
v = H * qd;  % spring velocity
T = H.' * (damp * v);  % torque on ankle due to damper
B = jacobian(T, qd);  % ankle damping due to damper
simplify(B - B.')  % equals zero for any parameter value

% TODO Derive the jacobian of K using H and HH. The following did not work

% K2 = stiff * H.' * (eye(3) - eabs*(eye(3)/sqrt(d.' * d) + d*(2*d.'*H)/(d.'*d)^(-3/2))) * H;
% K2 = -stiff*(H.'*H) - eabs*stiff*H.'*(d.'*d*eye(3) - d*d.')*H/((d.'*d)^(3/2));
% dabs = sqrt(d.' * d);
% K2 = -stiff*H.'*(1 + eabs/dabs - eabs*(d*d.')*H/dabs^3);
% K2 = -stiff*[HH(:,:,1).' * d, HH(:,:,2).' * d, HH(:,:,3).' * d] - stiff*(H.'*H);
% K2 = -stiff*[HH(:,:,1).' * (d-e), HH(:,:,2).' * (d-e), HH(:,:,3).' * (d-e)] - ...
%     stiff*H.'*(eye(3) + 

matlabFunction(d, H, K, B, 'File', 'mockup_ankle', 'vars', {q, a, b, stiff, damp, eabs})

%% Define an ankle mockup configuration and viz
% a = [[28e-3;65e-3;0], [-28e-3;65e-3;0], [0;65e-3;28e-3], [0;65e-3;-28e-3]];
% b = [[40e-3;-14e-3;0], [-40e-3;-14e-3;0], [0;-14e-3;41e-3], [0;-14e-3;-41e-3]];
% q = [0;0;0];
% stiff = [347,347,347,347];
% damp = [0,0,0,0];
% eabs = [1,1,1,1] * (70e-3 + 36e-3);  % natural length
% % compression =============
% by = -5e-3;
% br = 40e-3;
% ay = 10e-3;
% ar = 40e-3;
% a = [[ar;ay;0], [-ar;ay;0], [0;ay;ar], [0;ay;-ar]];
% b = [[br;by;0], [-br;by;0], [0;by;br], [0;by;-br]];
% q = [0;0;0];
% stiff = 300*[347,347,347,347];
% damp = [0,0,0,0];
% eabs = [1,1,1,1] * (20e-3);  % natural length

ay = 60e-3;
ar = 15e-3;
by = 10e-3; %-20e-3;
br = 55e-3;

a = [[ar;ay;0], [-ar;ay;0], [0;ay;ar], [0;ay;-ar]];
b = [[br;by;0], [-br;by;0], [0;by;br], [0;by;-br]];
q = [0;0;0];
stiff = 1*[347,347,347,347];
damp = [0,0,0,0];
eabs = [1,1,1,1] * (145e-3);  % natural length



Kankle = zeros(3);  % summed effect of all springs/dampers
Bankle = zeros(3);
for i = 1 : 4  % lopp tru each spring
    [~, ~, K, B] = mockup_ankle(q, a(:,i), b(:,i), stiff(i), damp(i), eabs(i));
    Kankle = Kankle + K;
    Bankle = Bankle + B;
end
Kankle, %Bankle

figure, plot_3d_ankle(q, a, b, stiff, damp, eabs);


