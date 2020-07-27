function hh = plot_3d_ankle(q, a, b, stiff, damp, eabs, hh)
% Input ==============================
%   x, y, q: end-effector pose
% Geometry ===========================
%   p11, p12, p13: inertial origin to passive rev joint on base of pistons [I frame]
%   beta1, beta2, beta3: orientation of pistons [I frame]
%   L1, L2, L3: length of extension links
%   r1, r2, r3: plate origin to passive rev joint [P frame]
% Calculated =========================
%   th1, th2, th3: orientation of extension links
%   p21, p22, p23: inertial origin to passive rev joint on end of pistons [I frame]
%   p31, p32, p33: inertial origin to passive rev joint on top plate [I frame]

% Assumptions ========================
%   - P and I frame are aligned when d1=d2=d3=0

arrow_len = 0.1;
num_z = size(a, 2);

% shank
p0 = [0; 0.3; 0];  % Useful later for the moving base
R0 = eye(3);

p = [0; 0; 0];
R = R0 * Ry(q(2)) * Rz(q(3)) * Rx(q(1));  % foot

pb = zeros(3, num_z); % b points in global frame
pba = zeros(3, num_z); % natural length
K = zeros(3);
B = zeros(3);
for i = 1 : num_z
    [d, H, Ki, Bi] = mockup_ankle(q, a(:,i), b(:,i), stiff(i), damp(i), eabs(i));
%     pb(:,i) = a(:,i) + Rz(Bet(i))*Rx(Alf(i)) * [0; D(i); 0];
%     pba(:,i) = a(:,i) + Rz(Bet(i))*Rx(Alf(i)) * [0; piston_len(i); 0];
    
    pb(:,i) = a(:,i) + R0 * d;
    pba(:,i) = pb(:,i) - R0 * eabs(i)*d/sqrt(d.'*d);
    K = K + Ki;
    B = B + Bi;
end
% K, B

feet = [0, 0, 0
        0, -3e-2, 0
        20e-2, -3e-2, 0] * R' + p';  % TODO Use patch instead
nperm = reshape(nchoosek(1:(num_z+1), 2)', [], 1);

% Correct =====================================
Rfix = [0, 1, 0
        0, 0, 1
        1, 0, 0]';
    
p0 = Rfix * p0;
p = Rfix * p;
pb = Rfix * pb;
pba = Rfix * pba;
a = Rfix * a;
feet = feet * Rfix';
R0 = Rfix * R0;
R = Rfix * R;

% plot ========================================

aa = [a, [0;0;0]];
pbb = [pb, p];

if nargin == 15
%     hh = [hpa, hW, hpb, hB, hleg, hmot, hfoot];
    hpa = hh(1); hW = hh(2:6); hpb = hh(7); hB = hh(8:12); hleg = hh(13:18);
    hmot = hh(19:24); hfoot = hh(25);
    
    
    set(hpa, 'xData', aa(1,nperm), 'yData', aa(2,nperm), 'zData', aa(3,nperm));
    plot_frame(p0, R0 * arrow_len, 'W', hW);
    set(hpb, 'xData', pbb(1,nperm), 'yData', pbb(2,nperm), 'zData', pbb(3,nperm));
    plot_frame(p, R * arrow_len, 'W', hB);
    
    set(hfoot, 'xData', feet(:,1), 'yData', feet(:,2), 'zData', feet(:,3));
    
    for i = 1 : num_z
        set(hleg(i), 'xData', [a(1,i), pb(1,i)], 'yData', [a(2,i), pb(2,i)], 'zData', [a(3,i), pb(3,i)]);
        set(hmot(i), 'xData', [pb(1,i), pba(1,i)], 'yData', [pb(2,i), pba(2,i)], 'zData', [pb(3,i), pba(3,i)]);
    end
    
    
%     set(h6, 'xData', x, 'yData', y, 'uData', [1,0] * Rz(q)*[arrow_len;0], 'vData', [0,1] * Rz(q)*[arrow_len;0]);
%     set(h7, 'xData', x, 'yData', y, 'uData', [1,0] * Rz(q)*[0; arrow_len], 'vData', [0,1] * Rz(q)*[0;arrow_len]);
%     set(h8, 'Position', [x+0.05, y+0.05]);
%     set(h8b, 'xData', x, 'YData', y);
%     
%     set(h15, 'xData', [p21a(1), p21(1)], 'YData', [p21a(2), p21(2)]);
% %     set(h16, 'position', [(p11+p21)/2; 0]', 'string', sprintf('  %.2f', d1-piston_len(1)));
%     set(h17, 'xData', [p22a(1), p22(1)], 'YData', [p22a(2), p22(2)]);
% %     set(h18, 'position', [(p12+p22)/2; 0]', 'string', sprintf('  %.2f', d2-piston_len(2)));
%     set(h19, 'xData', [p23a(1), p23(1)], 'YData', [p23a(2), p23(2)]);
% %     set(h20, 'position', [(p13+p23)/2; 0]', 'string', sprintf('  %.2f', d3-piston_len(3)));
%     
%     set(h21, 'xData', [p11(1), p21a(1)], 'YData', [p11(2), p21a(2)]);
%     set(h22, 'xData', [p12(1), p22a(1)], 'YData', [p12(2), p22a(2)]);
%     set(h23, 'xData', [p13(1), p23a(1)], 'YData', [p13(2), p23a(2)]);
%     
%     set(h24, 'xData', [p11(1), p21(1), p12(1), p22(1), p13(1), p23(1)], ...
%         'yData', [p11(2), p21(2), p12(2), p22(2), p13(2), p23(2)]);
%     
%     set(h25, 'xData', feet(:,1), 'yData', feet(:,2));
    
else
    
    hpa = plot3(aa(1,nperm), aa(2,nperm), aa(3,nperm), '-o');
    hold on
    hW = plot_frame(p0, R0 * arrow_len, 'W');

    hpb = plot3(pbb(1,nperm), pbb(2,nperm), pbb(3,nperm), '-o');
    hB = plot_frame(p, R * arrow_len, 'B');
    
    for i = 1 : num_z
        hleg(i) = plot3([a(1,i), pb(1,i)], [a(2,i), pb(2,i)], [a(3,i), pb(3,i)], 'k--');
        hmot(i) = plot3([pb(1,i), pba(1,i)], [pb(2,i), pba(2,i)], [pb(3,i), pba(3,i)], 'Color', [0,0,0, 0.5], 'linewidth', 5);
    end
    
    hfoot = plot3(feet(:,1), feet(:,2), feet(:,3), 'Color', [0.8235    0.6314    0.5490, 0.4], 'linewidth', 15);

    wbase = 0.3;
    xlim([-wbase, wbase]), ylim([-wbase, wbase]), zlim([-0.1, 0.4]) %grid on
%     xticks(-0.4:0.2:0.4), yticks(-0.0:0.2:0.6)
    set(gca, 'dataaspectratio', [1,1,1]), xlabel('z [m]'), ylabel('x [m]'), zlabel('y [m]')
    hold off
    set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')
    grid on
    
    hh = [hpa, hW, hpb, hB, hleg, hmot, hfoot];
end
end

function h = plot_frame(p, R, name, h)
    if nargin == 3
        hx = quiver3(p(1), p(2), p(3), R(1,1), R(2,1), R(3,1), 'r');
        hy = quiver3(p(1), p(2), p(3), R(1,2), R(2,2), R(3,2), 'g');
        hz = quiver3(p(1), p(2), p(3), R(1,3), R(2,3), R(3,3), 'b');
        h0 = plot3(p(1), p(2), p(3), 'sk', 'MarkerFaceColor', [0,0,0]);
        ht = text(p(1), p(2), p(3), sprintf('  %s', name));
        h = [hx, hy, hz, h0, ht];
    else
        set(h(1), 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'UData', R(1,1), 'VData', R(2,1), 'WData', R(3,1))
        set(h(2), 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'UData', R(1,2), 'VData', R(2,2), 'WData', R(3,2))
        set(h(3), 'XData', p(1), 'YData', p(2), 'ZData', p(3), 'UData', R(1,3), 'VData', R(2,3), 'WData', R(3,3))
        set(h(4), 'XData', p(1), 'YData', p(2), 'ZData', p(3))
        set(h(5), 'Position', p);
    end
end

    
    
    
    
    
