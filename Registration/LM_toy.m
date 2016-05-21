function [ F ] = LM_toy( P, data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


[n_unique, n_S, ~] = size(data);

   
    
S_R = cell(1,n_S);
S_t = cell(1,n_S);
ix = 1;
for s = 2:n_S;
    rx = P(ix);
    ry = P(ix+1);
    rz = P(ix+2);
    S_R{s} = compose_rotation(rx, ry, rz);
    S_t{s} = P(ix+3:ix+5);
    ix = ix + 6;
end
    

X_13_15 = [x(1:20);x(21:40);x(41:60)];
X_13_08 = [x(61:88);x(89:116);x(117:144)];
X_15_13 = [x(145:164);x(165:184);x(185:204)];
X_08_13 = [x(205:232);x(233:260);x(261:288)];
rx_15_13 = x(289);
ry_15_13 = x(290);
rz_15_13 = x(291);
rx_08_13 = x(292);
ry_08_13 = x(293);
rz_08_13 = x(294);
t_15_13 = x(295:297)';
t_08_13 = x(298:300)';

R_15_13 = compose_rotation(rx_15_13,ry_15_13,rz_15_13);
R_08_13 = compose_rotation(rx_08_13,ry_08_13,rz_08_13);

n_15 = size(X_15_13,2);
n_08 = size(X_08_13,2);
Y_15_13 = R_15_13 * X_15_13 + repmat(t_15_13,[1,n_15]);
Y_08_13 = R_08_13 * X_08_13 + repmat(t_08_13,[1,n_08]);

F = [(X_13_15 - Y_15_13) (X_13_08 - Y_08_13)];
F = reshape(F,[1,numel(F)]);

end