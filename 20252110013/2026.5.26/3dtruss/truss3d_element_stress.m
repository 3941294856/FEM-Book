function [epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de)
% 计算三维杆单元的轴向应变、应力和轴力
% 输入：
%   x1, x2 : 节点坐标 [x, y, z]
%   E, A   : 弹性模量 (Pa) 和截面积 (m^2)
%   de     : 单元节点位移列阵 [u1,v1,w1, u2,v2,w2] (m)
% 输出：
%   epsilon : 轴向应变 (无量纲)
%   sigma   : 轴向应力 (Pa)
%   N       : 轴力 (N)，拉正压负

    % 先计算长度和方向余弦
    [L, dir_cos, ~] = truss3d_element_stiffness(x1, x2, E, A);
    cx = dir_cos(1);
    cy = dir_cos(2);
    cz = dir_cos(3);
    
    % 节点位移
    u1 = de(1); v1 = de(2); w1 = de(3);
    u2 = de(4); v2 = de(5); w2 = de(6);
    
    % 轴向伸长量
    delta = (u2 - u1)*cx + (v2 - v1)*cy + (w2 - w1)*cz;
    
    % 应变
    epsilon = delta / L;
    
    % 应力 (胡克定律)
    sigma = E * epsilon;
    
    % 轴力
    N = sigma * A;
end