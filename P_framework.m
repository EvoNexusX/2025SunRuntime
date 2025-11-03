clear
clc
close all %清理橱窗
rng(1) %确保随机数生成是可重复的
addpath(genpath(pwd)); %将所有目录和子目录添加到运行文件夹下

PHI=0.95;
D_values = 10:10:200;

results_multi = zeros(1, length(D_values));
std_multi = zeros(1, length(D_values));
results_simple = zeros(1, length(D_values)); 
std_simple = zeros(1, length(D_values));
results_random = zeros(1, length(D_values));  % 预分配结果数组
std_random = zeros(1, length(D_values));
results_payoff = zeros(1, length(D_values)); 
std_payoff = zeros(1, length(D_values));
    
% 遍历每个D值，计算平均值
for i = 1:length(D_values)  %对于每一个问题大小
    dim = D_values(i);
    runs_multi= zeros(1, 10);
    runs_simple = zeros(1, 10);
    runs_random = zeros(1, 10); % 存储10次运行结果
    runs_payoff = zeros(1, 10);
    for j = 1:10
        runs_multi(j) = log(P_multiobjective(dim));
        runs_simple(j) = log(P_simple(dim));
        runs_random(j) = log(P_random_1(dim,PHI)); 
        runs_payoff(j) = log(P_payoff(dim));
    end
    results_multi(i) = mean(runs_multi); %对于每一个问题大小，得到10次运行的均值与方差
    std_multi(i) = std(runs_multi);
    results_simple(i) = mean(runs_simple);
    std_simple(i) = std(runs_simple);
    results_random(i) = mean(runs_random); % 计算10次运行的平均值
    std_random(i) = std(runs_random);
    results_payoff(i) = mean(runs_payoff);
    std_payoff(i) = std(runs_payoff);
end

% 绘制折线图
figure;
hold on; % 允许在同一图中绘制多条线
colors = lines(4); % 为不同问题生成颜色
legend_entries = {'multi','simple','random','payoff'}; % 图例内容

% 绘制每条线的浮动范围（使用fill绘制阴影区域）
plot(D_values, results_multi, '-', 'LineWidth', 1, 'Color', colors(1, :), 'DisplayName', 'multi');
plot(D_values, results_simple, '-', 'LineWidth', 1, 'Color', colors(2, :), 'DisplayName', 'simple');
plot(D_values, results_random, '-', 'LineWidth', 1, 'Color', colors(3, :), 'DisplayName', 'random');
plot(D_values, results_payoff, '-', 'LineWidth', 1, 'Color', colors(4, :), 'DisplayName', 'payoff');

fill([D_values, fliplr(D_values)], [results_multi + std_multi, fliplr(results_multi - std_multi)], colors(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
fill([D_values, fliplr(D_values)], [results_simple + std_simple, fliplr(results_simple - std_simple)], colors(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([D_values, fliplr(D_values)], [results_random + std_random, fliplr(results_random - std_random)], colors(3, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
fill([D_values, fliplr(D_values)], [results_payoff + std_payoff, fliplr(results_payoff - std_payoff)], colors(4, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 


% 添加网格、标签和标题
grid on;
xlabel('Problem size', 'FontSize', 12);
ylabel('The base 10 Average time', 'FontSize', 12);
% title('D值与平均结果的关系', 'FontSize', 14);

% 添加图例
legend(legend_entries, 'Location', 'Best');

% 释放 hold
hold off;

