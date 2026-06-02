%% 三维杆单元程序设计作业验证脚本
% 包含算例1、算例2以及刚度矩阵性质验证

clear; clc; format compact; format short g;

%% 算例1：沿x轴的一维杆单元
fprintf('========== 算例1：一维杆单元 ==========\n');
x1 = [0, 0, 0];
x2 = [2, 0, 0];
E = 200e9;      % 200 GPa
A = 1e-4;       % 1e-4 m^2
de = [0, 0, 0, 1e-3, 0, 0];  % 节点位移

[L, dc, Ke] = truss3d_element_stiffness(x1, x2, E, A);
[epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de);

fprintf('单元长度 L = %.6f m (理论值 2 m)\n', L);
fprintf('方向余弦 = (%.1f, %.1f, %.1f) 理论值 (1,0,0)\n', dc);
fprintf('刚度矩阵 Ke (6x6):\n');
disp(Ke);
fprintf('轴向应变 epsilon = %.4e (理论值 5e-4)\n', epsilon);
fprintf('轴向应力 sigma = %.2f MPa (理论值 100 MPa)\n', sigma/1e6);
fprintf('轴力 N = %.2f N (理论值 1e4 N)\n', N);

% 检查退化情况：位移只有x方向，应变应为 (u2-u1)/L = 1e-3/2 = 5e-4
fprintf('\n--- 验证完成 ---\n\n');

%% 算例2：空间任意方向杆单元
fprintf('========== 算例2：空间任意方向杆单元 ==========\n');
x1 = [0, 0, 0];
x2 = [1, 2, 2];
E = 210e9;      % 210 GPa
A = 2e-4;       % 2e-4 m^2
de = [0, 0, 0, 1e-3, 2e-3, 2e-3];

[L, dc, Ke] = truss3d_element_stiffness(x1, x2, E, A);
[epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de);

fprintf('单元长度 L = %.6f m (理论值 3 m)\n', L);
fprintf('方向余弦 = (%.4f, %.4f, %.4f) 理论值 (1/3, 2/3, 2/3)\n', dc);
fprintf('刚度矩阵 Ke (6x6):\n');
disp(Ke);

% 验证刚度矩阵性质
fprintf('\n--- 刚度矩阵性质验证 ---\n');
% 1. 对称性
if max(max(abs(Ke - Ke'))) < 1e-10
    fprintf('✓ 刚度矩阵对称\n');
else
    fprintf('✗ 刚度矩阵不对称\n');
end

% 2. 奇异性（行列式接近零）
det_Ke = det(Ke);
fprintf('行列式 det(Ke) = %.4e (应接近0)\n', det_Ke);

% 3. 半正定性（特征值非负）
eigvals = eig(Ke);
fprintf('特征值：\n');
for i = 1:length(eigvals)
    fprintf('  λ%d = %.4e\n', i, eigvals(i));
end
if all(eigvals >= -1e-10)
    fprintf('✓ 所有特征值非负，刚度矩阵半正定\n');
else
    fprintf('✗ 存在负特征值\n');
end
% 单个自由杆单元刚度矩阵奇异的原因：存在刚体位移模式（特征值=0）

% 4. 刚体位移验证：给单元施加一个刚体平移，应产生零内力
de_rigid = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3];  % 整体平移
[eps_r, sigma_r, N_r] = truss3d_element_stress(x1, x2, E, A, de_rigid);
fprintf('刚体平移位移 => 应变 = %.4e, 应力 = %.2f MPa, 轴力 = %.2f N (均应接近0)\n', ...
    eps_r, sigma_r/1e6, N_r);

% 计算应变、应力、轴力
fprintf('\n--- 单元应力计算结果 ---\n');
fprintf('轴向应变 epsilon = %.4e (理论值 1e-3)\n', epsilon);
fprintf('轴向应力 sigma = %.2f MPa (理论值 210 MPa)\n', sigma/1e6);
fprintf('轴力 N = %.2f N (理论值 4.2e4 N)\n', N);

%% 任务4：刚度矩阵物理意义验证
fprintf('\n========== 任务4：刚度矩阵物理意义验证 ==========\n');
% 选取自由度 j = 3（节点1的w方向位移），令 de = [0,0,1,0,0,0]
j = 3;
de_j = zeros(6,1);
de_j(j) = 1;
Fe = Ke * de_j;
fprintf('令第 %d 个自由度位移为1，其他为0，计算 Fe = Ke * de：\n', j);
fprintf('Fe = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]^T\n', Fe);
fprintf('该列向量正是 Ke 的第 %d 列：\n', j);
fprintf('Ke(:, %d) = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]^T\n', j, Ke(:,j));
fprintf('结论：k_{i,j} 的物理意义是：当第 j 个自由度产生单位位移而其他自由度位移为零时，在自由度 i 方向上产生的节点力。\n');

%% 附加题：简单空间桁架结构分析（可选）
fprintf('\n========== 附加题：简单空间桁架 ==========\n');
% 结构：四节点四面体桁架，节点1固定，节点2、3、4受集中力
% 节点坐标
nodes = [ 0,0,0;    % 节点1
          1,0,0;    % 节点2
          0,1,0;    % 节点3
          0,0,1 ];  % 节点4
% 单元连接 [单元编号, 节点i, 节点j]
elements = [1, 1, 2;
            2, 1, 3;
            3, 1, 4;
            4, 2, 3;
            5, 3, 4;
            6, 4, 2];
% 材料与截面
E_extra = 210e9;
A_extra = 1e-4;
% 载荷：节点2上施加 Fx = 1000 N，节点3上施加 Fy = -500 N，节点4上施加 Fz = 800 N
loads = [2, 1000, 0, 0; 3, 0, -500, 0; 4, 0, 0, 800];
% 约束：节点1完全固定
constraints = [1, 0, 0, 0];  % x,y,z 位移均为0

% 总自由度
n_node = size(nodes,1);
n_dof = 3;
NDOF = n_node * n_dof;
K_global = zeros(NDOF, NDOF);
F_global = zeros(NDOF, 1);

% 组装整体刚度矩阵
for i = 1:size(elements,1)
    n1 = elements(i,2);
    n2 = elements(i,3);
    x1_node = nodes(n1,:);
    x2_node = nodes(n2,:);
    [~, ~, Ke_local] = truss3d_element_stiffness(x1_node, x2_node, E_extra, A_extra);
    % 自由度映射
    dof = [ (n1-1)*3+1, (n1-1)*3+2, (n1-1)*3+3, ...
            (n2-1)*3+1, (n2-1)*3+2, (n2-1)*3+3 ];
    K_global(dof, dof) = K_global(dof, dof) + Ke_local;
end

% 施加载荷
for i = 1:size(loads,1)
    node_id = loads(i,1);
    F_global((node_id-1)*3+1) = loads(i,2);
    F_global((node_id-1)*3+2) = loads(i,3);
    F_global((node_id-1)*3+3) = loads(i,4);
end

% 处理约束（置大数法）
constrained_dof = [];
for i = 1:size(constraints,1)
    node_id = constraints(i,1);
    for d = 1:3
        if constraints(i, d+1) == 0
            constrained_dof = [constrained_dof, (node_id-1)*3 + d];
        end
    end
end
a = 1e12;
K_modified = K_global;
F_modified = F_global;
for idx = 1:length(constrained_dof)
    dof = constrained_dof(idx);
    K_modified(dof, :) = 0;
    K_modified(:, dof) = 0;
    K_modified(dof, dof) = a;
    F_modified(dof) = a * 0;  % 强制位移为0
end

% 求解位移
U_global = K_modified \ F_modified;

% 输出节点位移
fprintf('节点位移 (m):\n');
for i = 1:n_node
    ux = U_global((i-1)*3+1);
    uy = U_global((i-1)*3+2);
    uz = U_global((i-1)*3+3);
    fprintf('节点%d: (%.6e, %.6e, %.6e)\n', i, ux, uy, uz);
end

% 计算各杆应力与轴力
fprintf('\n单元应力与轴力:\n');
for i = 1:size(elements,1)
    n1 = elements(i,2);
    n2 = elements(i,3);
    x1_node = nodes(n1,:);
    x2_node = nodes(n2,:);
    de_local = [U_global((n1-1)*3+1), U_global((n1-1)*3+2), U_global((n1-1)*3+3), ...
                U_global((n2-1)*3+1), U_global((n2-1)*3+2), U_global((n2-1)*3+3)];
    [eps, sig, N_elem] = truss3d_element_stress(x1_node, x2_node, E_extra, A_extra, de_local);
    fprintf('单元%d: 应变=%.4e, 应力=%.2f MPa, 轴力=%.2f N\n', i, eps, sig/1e6, N_elem);
end