% ======================================================
% 有限元法示例：多边形逼近π
% ======================================================
clear; clc; close all;

%% 1. 基础参数与计算
n = [1, 2, 4, 8, 16, 32, 64, 128, 256];
h = 1 ./ n;
pi_exact = pi;
pi_n = n .* sin(pi ./ n);  % 基础近似值
e_n = abs(pi_exact - pi_n); % 基础误差

%% 2. 执行Wynn-ε外插计算
extrapolated_table = nan(size(pi_n));
extrapolated_table(3) = wynn_epsilon(pi_n(1:3));   % n=4（序列长度3）
extrapolated_table(5) = wynn_epsilon(pi_n(1:5));   % n=16（序列长度5）
extrapolated_table(7) = wynn_epsilon(pi_n(1:7));   % n=64（序列长度7）
extrapolated_table(9) = wynn_epsilon(pi_n(1:9));   % n=256（序列长度9）
e_extrapolated = abs(pi_exact - extrapolated_table); % 外插误差

%% 3. 输出表格
fprintf('=================================================================================\n');
fprintf('  n\t\t πₙ = n·sin(π/n)\t Wynn-ε外插结果\t 误差（外插后）\n');
fprintf('=================================================================================\n');
for i = 1:length(n)
    if isnan(extrapolated_table(i))
        fprintf('%4d\t %.15f\t %-15s\t %-15s\n', ...
            n(i), pi_n(i), '-', '-');
    else
        fprintf('%4d\t %.15f\t %.15f\t %.2e\n', ...
            n(i), pi_n(i), extrapolated_table(i), e_extrapolated(i));
    end
end
fprintf('精确π值（16位）：%.15f\n', pi_exact);
fprintf('=================================================================================\n');

%% 4. 绘制误差曲线（基础+外插）
figure('Color','w');
loglog(h, e_n, 'b-o', 'LineWidth',1.5, 'MarkerSize',6); hold on;
valid_idx = ~isnan(e_extrapolated);
loglog(h(valid_idx), e_extrapolated(valid_idx), 'r-^', 'LineWidth',1.5, 'MarkerSize',6);

% 拟合两条曲线的斜率
p1 = polyfit(log10(h), log10(e_n), 1); % 基础误差斜率（理论=2）
p2 = polyfit(log10(h(valid_idx)), log10(e_extrapolated(valid_idx)), 1); % 外插误差斜率

% 绘制拟合线
h_fit = logspace(log10(min(h)), log10(max(h)), 100);
loglog(h_fit, 10^(p1(2))*h_fit.^p1(1), 'b--', 'LineWidth',1.5);
loglog(h_fit, 10^(p2(2))*h_fit.^p2(1), 'r--', 'LineWidth',1.5);

% 图表标注
xlabel('h = 1/n（单元尺寸）');
ylabel('eₙ = |π - πₙ|（误差）');
title('多边形逼近π的收敛曲线（含Wynn-ε外插）');
grid on;
legend('基础近似误差', 'Wynn-ε外插误差', ...
    'Location','southeast');

fprintf('\n基础近似收敛阶：%.2f（理论值2）\n', p1(1));
fprintf('外插后收敛阶：%.2f\n', p2(1));

%% ==================== 函数定义 ====================
% Wynn-ε外插算法
function extrap_val = wynn_epsilon(sequence)
    N = length(sequence);
    epsilon = zeros(N+1, N+1); 
    % 初始化
    for n_idx = 0:N-1
        epsilon(1, n_idx+1) = 0;       
        epsilon(2, n_idx+1) = sequence(n_idx+1); 
    end
    % 递推计算
    for k = 1:N-1
        for n_idx = 0:N-1-k
            epsilon(k+2, n_idx+1) = epsilon(k, n_idx+2) + ...
                1/(epsilon(k+1, n_idx+2) - epsilon(k+1, n_idx+1));
        end
    end
    extrap_val = epsilon(end, 1); 
end