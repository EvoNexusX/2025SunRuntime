clear
clc
close all %清理橱窗
rng(1) %确保随机数生成是可重复的
addpath(genpath(pwd)); %将所有目录和子目录添加到运行文件夹下

maxFE_values = 1000:10:3000;

% 运行10次并计算均值和方差
num_runs = 10;
results_simple_1 = zeros(num_runs, length(maxFE_values));
results_simple_2 = zeros(num_runs, length(maxFE_values));
results_simple_3 = zeros(num_runs, length(maxFE_values));
results_ctoc_1 = zeros(num_runs, length(maxFE_values));
results_ctoc_2 = zeros(num_runs, length(maxFE_values));
results_ctoc_3 = zeros(num_runs, length(maxFE_values));
results_multi_1 = zeros(num_runs, length(maxFE_values));
results_multi_2 = zeros(num_runs, length(maxFE_values));
results_multi_3 = zeros(num_runs, length(maxFE_values));

for i = 1:num_runs
%     results_simple_1(i, :) = simple_1(maxFE_values);
%     results_simple_2(i, :) = simple_2(maxFE_values);
    results_simple_3(i, :) = simple_3(maxFE_values);
%     results_ctoc_1(i, :) = ctoc_1(maxFE_values);
%     results_ctoc_2(i, :) = ctoc_2(maxFE_values);
    results_ctoc_3(i, :) = ctoc_3(maxFE_values);
%     results_multi_1(i, :) = multi_1(maxFE_values);
%     results_multi_2(i, :) = multi_2(maxFE_values);
    results_multi_3(i, :) = multi_3(maxFE_values);
end

% 计算均值和方差
mean_simple_3 = mean(results_simple_3, 1);
std_simple_3 = std(results_simple_3, 0, 1);

mean_ctoc_3 = mean(results_ctoc_3, 1);
std_ctoc_3 = std(results_ctoc_3, 0, 1);

mean_multi_3 = mean(results_multi_3, 1);
std_multi_3 = std(results_multi_3, 0, 1);


figure;
hold on;
colors = lines(3); % 为不同问题生成颜色
legend_entries = {'simple','cons','DEMO'}; % 图例内容

% 画图
plot(maxFE_values, mean_simple_3, '-', 'Color', colors(1, :), 'LineWidth', 1);
plot(maxFE_values, mean_ctoc_3, '-', 'Color', colors(2, :), 'LineWidth', 1);
plot(maxFE_values, mean_multi_3, '-', 'Color', colors(3, :), 'LineWidth', 1);

% 绘制阴影（方差范围）
fill([maxFE_values, fliplr(maxFE_values)], ...
    [mean_simple_3 + std_simple_3, fliplr(mean_simple_3 - std_simple_3)], ...
    colors(1, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([maxFE_values, fliplr(maxFE_values)], ...
    [mean_ctoc_3 + std_ctoc_3, fliplr(mean_ctoc_3 - std_ctoc_3)], ...
    colors(2, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([maxFE_values, fliplr(maxFE_values)], ...
    [mean_multi_3 + std_multi_3, fliplr(mean_multi_3 - std_multi_3)], ...
    colors(3, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% 设置图形属性
xlabel('Fitness Evaluations', 'FontSize', 12);
ylabel('epsilon', 'FontSize', 12);
legend(legend_entries, 'Location', 'Best');
grid on;
hold off;


%% 简单算法

% 简单算法+UAV1
function results_simple_1 = simple_1(maxFE_values)
num_maxFE = length(maxFE_values);
results_simple_1 = zeros(1, num_maxFE); % 存储epsilon_simple结果
% simple
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_simple_1(1, i) = NP_simple0(maxFE, UAV1(1));
end
end

% 简单算法+UAV3
function results_simple_2 = simple_2(maxFE_values)
num_maxFE = length(maxFE_values);
results_simple_2 = zeros(1, num_maxFE); % 存储epsilon_simple结果
% simple
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_simple_2(1, i) = NP_simple0(maxFE, UAV3(1));
end
end

% 简单算法+UAV5
function results_simple_3 = simple_3(maxFE_values)
num_maxFE = length(maxFE_values);
results_simple_3 = zeros(1, num_maxFE); % 存储epsilon_simple结果
% simple
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_simple_3(1, i) = NP_simple0(maxFE, UAV5(1));
end
end


%% ctoc

% ctoc算法+UAV1
function results_ctoc_1 = ctoc_1(maxFE_values)
num_maxFE = length(maxFE_values);
results_ctoc_1 = zeros(1, num_maxFE); 
% ctoc
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_ctoc_1(1, i) = NP_ctoc(maxFE, UAV1(1));
end
end

% ctoc算法+UAV3
function results_ctoc_2 = ctoc_2(maxFE_values)
num_maxFE = length(maxFE_values);
results_ctoc_2 = zeros(1, num_maxFE); 
% ctoc
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_ctoc_2(1, i) = NP_ctoc(maxFE, UAV3(1));
end
end


% ctoc算法+UAV5
function results_ctoc_3 = ctoc_3(maxFE_values)
num_maxFE = length(maxFE_values);
results_ctoc_3 = zeros(1, num_maxFE); 
% ctoc
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_ctoc_3(1, i) = NP_ctoc(maxFE, UAV5(1));
end
end


%% 多目标算法

% multi算法+UAV1
function results_multi_1 = multi_1(maxFE_values)
num_maxFE = length(maxFE_values);
results_multi_1 = zeros(1, num_maxFE);
% multi
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_multi_1(1, i) = NP_multiobjective(maxFE, UAV1(1));
end
end

% multi算法+UAV3
function results_multi_2 = multi_2(maxFE_values)
num_maxFE = length(maxFE_values);
results_multi_2 = zeros(1, num_maxFE);
% multi
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_multi_2(1, i) = NP_multiobjective(maxFE, UAV3(1));
end
end

% multi算法+UAV5
function results_multi_3 = multi_3(maxFE_values)
num_maxFE = length(maxFE_values);
results_multi_3 = zeros(1, num_maxFE); 
% multi
for i = 1:num_maxFE
    maxFE = maxFE_values(i);
    results_multi_3(1, i) = NP_multiobjective(maxFE, UAV5(1));
end
end