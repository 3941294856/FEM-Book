function F = Load_cases(Loads,F,Node_Dof)
%   PPT上抄的
%   将输入的载荷信息矩阵（Load）中的载荷提取出来组装成载荷列向量F
[row,col]=size(Loads);
for i=1:row
    for j=2:col
        F(Node_Dof*(Loads(i,1)-1)+j-1,1)=Loads(i,j);
    end
end