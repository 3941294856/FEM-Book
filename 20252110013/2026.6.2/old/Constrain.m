% ********************************************************************** %
% 引入约束条件 %
% ********************************************************************** %
% 获取单元类型、约束列表以及单元节点对应的自由度数
function [Constr_Pos,Constr_Num] = Constrain(BCs, Node_Dof)
%若为平面杆单元，取约束列表中的前3列数据，若为平面梁单元，取约束列表中的前4列数据
Constrain = BCs(:,1:(Node_Dof+1)); 
[row, col] = size(Constrain); % 获取新的约束列表的行、列数
m = 1;
for i = 1:row % 读取第一行数据
 for j = 2:col % 对各列进行遍历
 if Constrain(i,j) == 0; % 找到‘0’元素所在的行、列位置
 Constr_Dof(m,1) = Node_Dof * (Constrain(i,1)-1)+j-1;
 m = m + 1;
 end
 end
end
Constr_Pos = Constr_Dof;
Constr_Num = length(Constr_Dof);