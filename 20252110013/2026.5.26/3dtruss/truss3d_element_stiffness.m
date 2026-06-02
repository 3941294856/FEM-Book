function [L, dir_cos, Ke] = truss3d_element_stiffness(x1, x2, E, A)
% 计算三维杆单元的长度、方向余弦和全局刚度矩阵
% 输入：
%   x1, x2 : 节点坐标，形如 [x, y, z]
%   E      : 弹性模量 (Pa)
%   A      : 截面积 (m^2)
% 输出：
%   L          : 单元长度 (m)
%   dir_cos    : 方向余弦 [cx, cy, cz]
%   Ke         : 6x6 单元刚度矩阵 (N/m)

    % 计算坐标差
    dx = x2(1) - x1(1);
    dy = x2(2) - x1(2);
    dz = x2(3) - x1(3);
    
    L = sqrt(dx^2 + dy^2 + dz^2);
    
    % 退化单元检查
    if L < 1e-12
        error('错误：两个节点重合，单元长度为零，无法继续计算。');
    end
    
    % 方向余弦
    cx = dx / L;
    cy = dy / L;
    cz = dz / L;
    dir_cos = [cx, cy, cz];
    
    % 计算刚度矩阵
    k = E * A / L;
    % 6x6 刚度矩阵 (对称)
    Ke = k * [ cx*cx, cx*cy, cx*cz, -cx*cx, -cx*cy, -cx*cz;
               cx*cy, cy*cy, cy*cz, -cx*cy, -cy*cy, -cy*cz;
               cx*cz, cy*cz, cz*cz, -cx*cz, -cy*cz, -cz*cz;
              -cx*cx, -cx*cy, -cx*cz,  cx*cx,  cx*cy,  cx*cz;
              -cx*cy, -cy*cy, -cy*cz,  cx*cy,  cy*cy,  cy*cz;
              -cx*cz, -cy*cz, -cz*cz,  cx*cz,  cy*cz,  cz*cz ];
end