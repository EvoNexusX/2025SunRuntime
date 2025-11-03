%检查random的理论分析
% 先找到一方的一个PF，再进化另外一方，去寻找公共解,这个时间量级应该比checksimple的要多一点点，加个1/2没区别，数量级是一样的，这个正好对应理论分析的过程
%他这个没办法种群中只剩公共解，因为不论你怎么筛，种群中总会有PF2中的个体和公共解
%随机的那个因为他的方向每一轮是变化的，然后通过删去种群中被支配的粒子可以只剩下公共解
%这个你要是加上删除种群中被支配的个体纯属多余，因为你这个是两个方向单独分开来的，没用还多余

%测试正确

% 
% clear
% clc
% close all %清理橱窗
% rng(1) %确保随机数生成是可重复的
% addpath(genpath(pwd)); %将所有目录和子目录添加到运行文件夹下

function t = P_checkrandom(D,PHI)
%% the hole framework begin
% D = 200; % 问题维度
dm = 2; % 参与者数量，一共有两个人
nobj = 4; %目标数量,一共有4个目标
singlenobj = 2; % 每个参与方的有两个目标
common_pf = ones(1,D); %公共解(全1向量)
common_obj = BPAOAZ(common_pf); %公共解对应的目标值
common = [common_pf common_obj];

%第一方的pf
first_col = (0:D/2)'; %从0到n/2
second_col = (D:-1:D/2)'; %从n/2到n   这里因为默认n都是偶数，所以不做取整操作
pf_1 = [first_col, second_col];


%% initialize the population P
X0 = zeros(1,D); %使用全0向量初始化种群
Y0 = BPAOAZ(X0); % objectives of X0

population = [X0 Y0]; % 初始种群
%% begin iterations to uodate population
tic; % 开始计时

while  ~any(ismember(population(:,D+1:D+singlenobj), pf_1,  'rows')) %先向第一方优化，找到PF1中的一个个体
    %随机均匀地从种群P中选择一个粒子
    randnum = randi([1, size(population,1)]); % 从1到种群粒子数量之间随机选择一个整数
    mutate_x = population(randnum,1:D); %从种群中随机选择到的一个粒子
    x = mutate_x;%待变异的粒子

    choosed_direction = 0; % 此时优化方向为参与方1
    
    %选择完方向之后，粒子进行进化 one-bit mutation，即随机选择一位然后将其变异
    mutation_index = randi([1, D]); % 随机选择一个位置
    x(mutation_index) = 1-x(mutation_index); %完成一位翻转
    obj_x = BPAOAZ(x); %计算粒子的目标值
    xs =[x  obj_x];

    dominated = false;
    for z = 1:size(population,1) %这里判断的时候用弱支配
        if weakdominates(population(z,D+1+choosed_direction:D+singlenobj+choosed_direction), obj_x(1+choosed_direction:singlenobj+choosed_direction))
            dominated = true;%这时就说明种群中存在比他好或者跟他一样的粒子，这时就不用添加新粒子了。
            break;
        end
    end

    if ~dominated %说明种群中没有比他好的，可以更新种群，去掉种群中被他支配的粒子
        population = [population(~arrayfun(@(z) dominates(... %有相等的粒子的话在上一步就被卡掉了。能到这一步说明没有相等的，用大于等于没事
            obj_x(1+choosed_direction:singlenobj+choosed_direction), population(z,D+1+choosed_direction:D+singlenobj+choosed_direction)), 1:size(population,1)), :); xs];
    end
    
end



while ~ismember(common, population, 'rows') %找到PF1中的一个个体之后再去向第二方进化，直到找到公共解停止
    %随机均匀地从种群P中选择一个粒子
    randnum = randi([1, size(population,1)]); % 从1到种群粒子数量之间随机选择一个整数
    mutate_x = population(randnum,1:D); %从种群中随机选择到的一个粒子
    x = mutate_x;%待变异的粒子

    choosed_direction = singlenobj; % 此时优化方向为参与方2
    
    %选择完方向之后，粒子进行进化 one-bit mutation，即随机选择一位然后将其变异
    mutation_index = randi([1, D]); % 随机选择一个位置
    x(mutation_index) = 1-x(mutation_index); %完成一位翻转
    obj_x = BPAOAZ(x); %计算粒子的目标值
    xs =[x  obj_x];

    dominated = false;
    for z = 1:size(population,1) %这里判断的时候用弱支配
        if weakdominates(population(z,D+1+choosed_direction:D+singlenobj+choosed_direction), obj_x(1+choosed_direction:singlenobj+choosed_direction))
            dominated = true;%这时就说明种群中存在比他好或者跟他一样的粒子，这时就不用添加新粒子了。
            break;
        end
    end

    if ~dominated %说明种群中没有比他好的，可以更新种群，去掉种群中被他支配的粒子
        population = [population(~arrayfun(@(z) dominates(... %有相等的粒子的话在上一步就被卡掉了。能到这一步说明没有相等的，用大于等于没事
            obj_x(1+choosed_direction:singlenobj+choosed_direction), population(z,D+1+choosed_direction:D+singlenobj+choosed_direction)), 1:size(population,1)), :); xs];
    end

%     if size(population,1) > 1
%         population = removeNondominated(population,choosed_direction,singlenobj,D);
%     end
    
end
t = toc; % 结束计时并返回运行时间
% fprintf('random算法dim=%.0f时运行时间为 %.7f 秒\n', D,t);

end




%% 支配关系的定义
function is_dominated = dominates(x, y) %支配，不包含相等的情况，至少有一个不同
    is_dominated = false;
    if all(x >= y) && any(x > y) 
        is_dominated = true;
    end
end

function weak_dominated = weakdominates(x, y) %弱支配的情况，包含相等
    weak_dominated = false;
    if all(x >= y)
        weak_dominated = true;
    end
end

% function newpopulation = removeNondominated(population,d,singlenobj,D) %删除种群中在这个方向上的非支配解
%     n = size(population, 1); % 种群大小
%     is_dominated = false(n, 1); % 标记每个解是否被支配
% 
%     % 遍历种群中的每一对粒子
%     for i = 1:n
%         for j = 1:n %对于第i个粒子，判断他是不是被种群中的其他粒子支配
%             if i ~= j 
%                 % 如果粒子j支配粒子i，则标记i为被支配的粒子
%                 if weakdominates(population(j,D+1+d:D+singlenobj+d), population(i,D+1+d:D+singlenobj+d))
%                     is_dominated(i) = true;
%                     break; % 一旦被支配，无需再检查
%                 end
%             end
%         end
%     end
% 
%     % 过滤出非支配解
%     newpopulation = population(~is_dominated, :);
% 
% end



%% 两方两目标问题
function result = BPAOAZ(x)
    % BPAOAZ: Computes the pseudo-Boolean function BPAOAZ.
    % Input:
    %   x - A binary vector of length n, where n is even.
    % Output:
    %   AORZ - A 2-element vector [f_11(x), f_12(x)].
    %   AOFZ - A 2-element vector [f_21(x), f_22(x)].
  
    n = length(x);
    % Split the vector into two halves.
    n_half = n / 2;
    x1 = x(1:n_half);        % First half of x.
    x2 = x(n_half+1:end);    % Second half of x.
    
    % Compute f_11 and f_12 for AORZ.
    f_11 = sum(x2); %后1
    f_12 = sum(x1) + sum(1 - x2); %前1+后0
    AORZ = [f_11, f_12];
    
    % Compute f_21 and f_22 for AOFZ.
    f_21 = sum(1 - x1) + sum(x2); %前0+后1
    f_22 = sum(x1); %前1
    AOFZ = [f_21, f_22];
    
    result =[f_11, f_12, f_21, f_22];
end


%% 检查种群是否满足要求
function isEqual = check_1(P, D, singlenobj, pf) %pf_1
    sortedP = sortrows(P, D+1);
    isEqual = isequal(sortedP(:,D+1:D+singlenobj), pf);
end

function isEqual = check_2(P, D, singlenobj, pf) %pf_2
    sortedP = sortrows(P, D+singlenobj+singlenobj);
    isEqual = isequal(sortedP(:,D+1+singlenobj:D+singlenobj+singlenobj), pf);
end
