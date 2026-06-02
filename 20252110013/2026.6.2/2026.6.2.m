% **********************************************************************
% 程序设计作业: 二维桁架结构总体刚度矩阵组装与求解
% 符合《2.3 系统总体刚度方程》作业要求
% 功能: 前处理、单元分析、LM矩阵生成、直接组装、缩减法求解、后处理
% **********************************************************************

clear; clc; close all;

%% ====================== 选择算例 ======================
% 1: 一维两单元杆结构(算例1)
% 2: 二维两杆桁架结构(算例2)
example = 2;

%% ====================== 1. 前处理模块 ======================
fprintf('==================== 前处理 ====================\n');

if example == 1
    % 算例1: 一维两单元杆结构
    Title = '一维两单元杆结构';
    nsd = 1;          % 空间维数
    ndof = 1;         % 每个节点自由度
    nnp = 3;          % 节点总数
    nel = 2;          % 单元总数
    nen = 2;          % 每个单元节点数
    
    % 节点坐标(列: x)
    Nodes = [0; 1; 2];
    
    % 单元连接数组IEN(行: 单元, 列: 节点编号, 从1开始)
    IEN = [1, 2; 2, 3];
    
    % 材料参数
    E = [100*1, 200*1];  % E*A/L 直接给出
    A = [1, 1];          % 截面积(此处为1, 因为E*A/L已给出)
    
    % 边界条件: 已知位移自由度(从1开始)和位移值
    fixed_dof = [1];
    fixed_value = [0.0];
    
    % 节点载荷: 载荷自由度(从1开始)和载荷值
    force_dof = [2, 3];
    force_value = [0.0, 10.0];
    
elseif example == 2
    % 算例2: 二维两杆桁架结构
    Title = '二维两杆桁架结构';
    nsd = 2;          % 空间维数
    ndof = 2;         % 每个节点自由度
    nnp = 3;          % 节点总数
    nel = 2;          % 单元总数
    nen = 2;          % 每个单元节点数
    
    % 节点坐标(列: x, y)
    Nodes = [1.0, 0.0;  % 节点1
             0.0, 0.0;  % 节点2
             1.0, 1.0]; % 节点3
    
    % 单元连接数组IEN(行: 单元, 列: 节点编号, 从1开始)
    IEN = [1, 3; 2, 3];
    
    % 材料参数
    E = [1.0, 1.0];    % 弹性模量
    A = [1.0, 1.0];    % 截面积
    
    % 边界条件: 已知位移自由度(从1开始)和位移值
    % 自由度顺序: [u1, v1, u2, v2, u3, v3]
    fixed_dof = [1, 2, 3, 4];
    fixed_value = [0.0, 0.0, 0.0, 0.0];
    
    % 节点载荷: 载荷自由度(从1开始)和载荷值
    force_dof = [5, 6];
    force_value = [10.0, 0.0];
end

% 总自由度
DOFs = ndof * nnp;
fprintf('算例: %s\n', Title);
fprintf('节点数: %d, 单元数: %d, 总自由度: %d\n', nnp, nel, DOFs);

%% ====================== 2. 生成对号矩阵LM ======================
fprintf('\n==================== 生成对号矩阵LM ====================\n');

% LM矩阵: 行=单元自由度, 列=单元, 值=总体自由度编号(从1开始)
LM = zeros(nen*ndof, nel);

for e = 1:nel
    for i = 1:nen
        node_id = IEN(e, i);  % 单元第i个节点的全局编号
        for d = 1:ndof
            % 计算单元自由度在LM中的行号
            lm_row = (i-1)*ndof + d;
            % 计算对应的总体自由度编号
            LM(lm_row, e) = (node_id-1)*ndof + d;
        end
    end
end

fprintf('对号矩阵LM:\n');
disp(LM);

%% ====================== 3. 单元分析与总体刚度矩阵组装 ======================
fprintf('\n==================== 总体刚度矩阵组装 ====================\n');

% 初始化总体刚度矩阵
K = zeros(DOFs, DOFs);

% 存储每个单元的信息
elem_info = struct('length', [], 'c', [], 's', [], 'Ke', []);

for e = 1:nel
    % 获取单元节点坐标
    node1 = IEN(e, 1);
    node2 = IEN(e, 2);
    x1 = Nodes(node1, 1);
    y1 = Nodes(node1, 2);
    x2 = Nodes(node2, 1);
    y2 = Nodes(node2, 2);
    
    % 计算单元长度和方向余弦
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    c = (x2 - x1) / L;
    s = (y2 - y1) / L;
    
    % 计算单元刚度矩阵
    if nsd == 1
        % 一维杆单元
        Ke = E(e)*A(e)/L * [1, -1; -1, 1];
    else
        % 二维桁架单元
        Ke = E(e)*A(e)/L * [c^2, c*s, -c^2, -c*s;
                            c*s, s^2, -c*s, -s^2;
                           -c^2, -c*s, c^2, c*s;
                           -c*s, -s^2, c*s, s^2];
    end
    
    % 存储单元信息
    elem_info(e).length = L;
    elem_info(e).c = c;
    elem_info(e).s = s;
    elem_info(e).Ke = Ke;
    
    % 直接组装算法: K[LM[a,e], LM[b,e]] += Ke[a,b]
    for a = 1:size(Ke,1)
        for b = 1:size(Ke,2)
            K(LM(a,e), LM(b,e)) = K(LM(a,e), LM(b,e)) + Ke(a,b);
        end
    end
end

fprintf('总体刚度矩阵K:\n');
disp(K);

% 检查对称性
is_symmetric = isequal(K, K');
fprintf('总体刚度矩阵是否对称: %s\n', mat2str(is_symmetric));

% 检查奇异性(施加边界条件前)
rank_before = rank(K);
fprintf('施加边界条件前K的秩: %d, 总自由度: %d\n', rank_before, DOFs);
fprintf('施加边界条件前K是否奇异: %s\n', mat2str(rank_before < DOFs));

%% ====================== 4. 边界条件处理与方程求解(缩减法) ======================
fprintf('\n==================== 边界条件处理与方程求解 ====================\n');

% 初始化载荷向量和位移向量(全部为列向量)
F = zeros(DOFs, 1);
U = zeros(DOFs, 1);

% 施加节点载荷
for i = 1:length(force_dof)
    F(force_dof(i)) = force_value(i);
end

% 划分已知位移自由度(E)和未知位移自由度(F)
E_dof = fixed_dof(:);          % 强制转为列向量
F_dof = setdiff(1:DOFs, E_dof);  % 未知位移自由度(列向量)

U_E = fixed_value(:);  % 关键修复: 强制转为列向量(n×1)

% 划分刚度矩阵子块
K_FF = K(F_dof, F_dof);
K_EF = K(E_dof, F_dof);
K_FE = K_EF';
K_EE = K(E_dof, E_dof);

F_F = F(F_dof);  % 未知位移对应的载荷(列向量)

% 求解未知位移: K_FF * U_F = F_F - K_FE * U_E
% 现在维度匹配: K_FE(m×n) * U_E(n×1) = m×1, 与F_F(m×1)相减
U_F = K_FF \ (F_F - K_FE * U_E);

% 重构完整位移向量
U(F_dof) = U_F;
U(E_dof) = U_E;

% 计算约束反力: R_E = K_EE * U_E + K_EF * U_F - F_E
F_E = F(E_dof);
R_E = K_EE * U_E + K_EF * U_F - F_E;

fprintf('节点位移U:\n');
disp(U);

fprintf('约束反力R_E(对应自由度%d):\n', E_dof);
disp(R_E);

% 检查缩减矩阵是否非奇异
rank_after = rank(K_FF);
fprintf('缩减矩阵K_FF的秩: %d, 未知自由度数: %d\n', rank_after, length(F_dof));
fprintf('缩减矩阵K_FF是否非奇异: %s\n', mat2str(rank_after == length(F_dof)));

%% ====================== 5. 后处理: 单元应力和轴力计算 ======================
fprintf('\n==================== 单元应力和轴力计算 ====================\n');

for e = 1:nel
    % 提取单元节点位移
    de = U(LM(:,e));
    
    L = elem_info(e).length;
    c = elem_info(e).c;
    s = elem_info(e).s;
    
    % 计算单元应力
    if nsd == 1
        % 一维杆单元: sigma = E/L * [-1, 1] * de
        sigma = E(e)/L * [-1, 1] * de;
    else
        % 二维桁架单元: sigma = E/L * [-c, -s, c, s] * de
        sigma = E(e)/L * [-c, -s, c, s] * de;
    end
    
    % 计算单元轴力
    N = sigma * A(e);
    
    % 存储结果
    elem_info(e).sigma = sigma;
    elem_info(e).N = N;
    
    % 输出单元信息
    fprintf('单元%d:\n', e);
    fprintf('  长度: %.6f\n', L);
    if nsd == 2
        fprintf('  方向余弦: c=%.6f, s=%.6f\n', c, s);
    end
    fprintf('  应力: %.6f\n', sigma);
    fprintf('  轴力: %.6f\n', N);
end

%% ====================== 6. 结果验证 ======================
fprintf('\n==================== 结果验证 ====================\n');

if example == 1
    % 算例1验证
    expected_K = [100, -100, 0; -100, 300, -200; 0, -200, 200];
    expected_U = [0; 0.1; 0.15];
    expected_R = 10;
    
    K_error = max(max(abs(K - expected_K)));
    U_error = max(abs(U - expected_U));
    R_error = abs(R_E - expected_R);
    
    fprintf('总体刚度矩阵最大误差: %.10f\n', K_error);
    fprintf('节点位移最大误差: %.10f\n', U_error);
    fprintf('约束反力误差: %.10f\n', R_error);
    
elseif example == 2
    % 算例2验证
    expected_u3 = 38.284271;
    expected_v3 = -10.000000;
    expected_sigma1 = -10.000000;
    expected_sigma2 = 14.142136;
    
    u3_error = abs(U(5) - expected_u3);
    v3_error = abs(U(6) - expected_v3);
    sigma1_error = abs(elem_info(1).sigma - expected_sigma1);
    sigma2_error = abs(elem_info(2).sigma - expected_sigma2);
    
    fprintf('节点3 u位移误差: %.10f\n', u3_error);
    fprintf('节点3 v位移误差: %.10f\n', v3_error);
    fprintf('单元1应力误差: %.10f\n', sigma1_error);
    fprintf('单元2应力误差: %.10f\n', sigma2_error);
end

fprintf('\n程序运行完成!\n');