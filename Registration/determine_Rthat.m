function [ Ptrue, Phat] = determine_Rthat( P_LCSi, P_radi, d_rx, d_ry, d_rz, d_tx, d_ty, d_tz, action )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_pts = size(P_LCSi,2);
P_radj = P_radi;
% Adjust rx
n_i = numel(d_rx);
Phat = zeros(n_i,6);
if ~action.rx;
    d_rx = 0*d_rx;
end
if ~action.ry;
    d_ry = 0*d_ry;
end
if ~action.rz;
    d_rz = 0*d_rz;
end
if ~action.tx;
    d_tx = 0*d_tx;
end
if ~action.ty;
    d_ty = 0*d_ty;
end
if ~action.tz;
    d_tz = 0*d_tz;
end

Ptrue = [d_rx; d_ry; d_rz; d_tx; d_ty; d_tz];

for i = 1:n_i;
    % Initial truth transformation parameters
    rx = d_rx(i);
    ry = d_ry(i);
    rz = d_rz(i);
    tx = d_tx(i);
    ty = d_ty(i);
    tz = d_tz(i);
    % Transform P_LCSi into P_LCSj
    R = compose_rotation(deg2rad(rx),deg2rad(ry),deg2rad(rz));
    t = [tx ty tz]';
    P_LCSj = R*P_LCSi + repmat(t,[1,n_pts]);
    [Rhat,that, ~, ~ ] = toy_registrationfunction(P_LCSi',P_LCSj',P_radi,P_radj);
    Phat(1,i) = rad2deg(decompose_rotation_rx(Rhat));
    Phat(2,i) = rad2deg(decompose_rotation_ry(Rhat));
    Phat(3,i) = rad2deg(decompose_rotation_rz(Rhat));
    Phat(4,i) = that(1);
    Phat(5,i) = that(2);
    Phat(6,i) = that(3);
end

end

