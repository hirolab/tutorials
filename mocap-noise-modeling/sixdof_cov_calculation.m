%% Calculate covariance matrix of mocap
% Add mocap-utils to path

%% Generate data

% num_markers = 30;
% r = 0.1 * randn(num_markers, 3);

r = [0, 0.2, 0.07
    0, 0.2, -0.07
    0, -0.2, 0.04
    0, -0.2, -0.04
    0.05, 0.15, 0
    0.05, -0.15, 0];  % markers on body frame. Shank
% r = [0.25, -0.1, 0.07
%     0.25, -0.1, -0.07
%     -0.1, -0.05, 0
%     0.15, -0.05, 0];  % markers on body frame. Foot
% r = [0.25, 0, 0.15
%     0.25, 0, -0.15
%     -0.25, 0, -0.15
%     -0.25, 0, 0.15
%     0.25, 0, -0.03
%     -0.25, 0, -0.03
%     0.15, 0, 0.15
%     -0.05, 0, 0.15
%     0, 0, -0.15];  % markers on body frame. Plate

num_markers = size(r, 1);
num_samples = 10000;

trans = [0, 0, 0] .* ones(num_samples, 3);  % nominal body position
quat = angle2quat(0, 0, 0, 'YZX') .* ones(num_samples, 4);  % nominal body angle

% create 6dof struct compatible with QTMParser
marker = repmat(struct('trans', 0, 'name', ''), [num_markers, 1]);
for i = 1 : num_markers
    marker(i).trans = trans + quatrotate(quatinv(quat), r(i,:));
    marker(i).name = sprintf('marker %d', i);
end
b = struct('name', 'body', 'trans', trans, 'quat', quat, 'marker', marker);

figure, QTMParser.plotBody(b)

% Validate transformation
% [transf, quatf] = QTMParser.recompute_6dof(b);
% figure, plot([trans-transf])
% figure, plot([quat-quatf])

% print -dpng imgs/plate.png -r200

%% Add noise in each marker and channel

markerTransCov = 0.1e-3 ^ 2;

bnoise = b;  % body with noise
for i = 1 : num_markers
    bnoise.marker(i).trans = b.marker(i).trans + ...
        sqrt(markerTransCov) * randn(num_samples, 3);
end
[transf, quatf] = QTMParser.recompute_6dof(bnoise);
qf = log_quat(quatmultiply(quatinv(quat), quatf));

qCov = cov(qf)
transCov = cov(trans - transf)
poseCov = cov([qf * 180/pi, (trans - transf) * 1e3])  % in deg and mm for similar units

fprintf('95%% Confidence Intervals: \n Angle: [%.3f, %.3f, %.3f]^T deg\n Position: [%.3f, %.3f, %.3f]^T mm\n\n', ...
    1.96*sqrt(diag(qCov))*180/pi, 1.96*sqrt(diag(transCov))*1e3)

figure, set(gcf, 'position', [680   733   433   365])
imagesc(abs(poseCov)), colormap('gray'), colorbar, caxis([0, 10e-3])
title('Absolute Pose Covariance [deg] \times [mm]'), xticklabels({'q_x', 'q_y', 'q_z', 't_x', 't_y', 't_z'})
yticklabels({'q_x', 'q_y', 'q_z', 't_x', 't_y', 't_z'})
axis square

% print -dpng imgs/shankCov.png -r200
