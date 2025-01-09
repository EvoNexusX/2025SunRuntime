clear
clc
close all %清理橱窗
rng(1) %确保随机数生成是可重复的
addpath(genpath(pwd)); %将所有目录和子目录添加到运行文件夹下

PHI_values=0.1:0.01:0.9;
D_values = 50:50:200;
num_runs = 10; % 每个phi的运行次数

results_random = zeros(length(D_values), length(PHI_values)); % 预分配结果数组
std_random = zeros(length(D_values), length(PHI_values));

for d = 1:length(D_values)
    D = D_values(d);
    for i = 1:length(PHI_values)
        PHI = PHI_values(i);
        runs_random = zeros(1, num_runs); % 存储num_runs运行结果
        for j = 1:num_runs
            runs_random(j) = P_random_1(D,PHI); 
        end
        results_random(d,i) = mean(runs_random); % 计算10次运行的平均值
        std_random(d,i) = std(runs_random);
    end
end


% 绘制折线图
figure;
hold on; % 允许在同一图中绘制多条线
colors = lines(4); 
legend_entries = {'dim=50','dim=100','dim=150','dim=200'}; % 图例内容
for d = 1:length(D_values)
    % 绘制当前维度的折线并添加图例
    plot(PHI_values, results_random(d, :),'-', 'Color', colors(d, :), 'LineWidth', 1);
end

for d = 1:length(D_values)
    fill([PHI_values, fliplr(PHI_values)], ...
         [results_random(d, :) + std_random(d, :), fliplr(results_random(d, :) - std_random(d, :))], colors(d, :), ...
         'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
end

% 添加网格、标签和标题
xlabel('Phi', 'FontSize', 12);
ylabel('Average time', 'FontSize', 12);

% 添加图例
legend(legend_entries, 'Location', 'Best');

grid on;
hold off;

