% *************************************************************************** %
% 该函数将引入边界条件，对平衡方程进行修正组 %
% 采用置大数法对整体刚度矩阵及载荷列向量进行修正 %
% 方法：在被约束的自由度，将其所对应的刚度矩阵对角线元素乘以一个较大的数 a %
% 将刚度矩阵元素K(i, i)变为aK(i, i),在对应的载荷处乘以系数（aK(i, i)） %
% **************************************************************************** %
function [Km,Fm] = modified_equation (K,F,Constr_Num,Constr_Dof)
a = 1.E10;
for i =1:Constr_Num
 K(Constr_Dof(i,1),Constr_Dof(i,1)) = a * K(Constr_Dof(i,1),Constr_Dof(i,1));
 F(Constr_Dof(i,1)) = K(Constr_Dof(i,1),Constr_Dof(i,1)) * F(Constr_Dof(i,1));
end
Km = K;
Fm = F;