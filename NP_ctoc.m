% clear
% clc
% close all
% rng(1)
% addpath(genpath(pwd));

function final_epsilon = NP_ctoc(maxFE,problem)
% 初始化种群
% maxFE = 1000;
% problem = UAV5(1);
epsilon = problem.epsilon; %近似程度，两方相同
D = problem.D;
r = epsilon^(1/D);

dm = 2; % 参与者数量
nobj = 4; %目标数量
singlenobj = 2;

ub = problem.points; %顶点的总数量
lb = 2;

X0 = [1, randi([lb, ub])];
Y0 = problem.CalObj(X0);


population = {X0};  
obj_population = {Y0}; 

%公共解以及公共解对应的目标值
common = problem.common;
obj_common = zeros(size(common,1),nobj);
for c = 1:size(common,1)
    obj_common(c,:) = problem.CalObj(common(c,:));
end


%开始变异
while problem.calcount <= maxFE
    %fprintf("FE:%d/MaxFE:%d\n",problem.calcount,maxFE);
    %随机均匀地从种群P中选择一个粒子
    randnum = randi([1, numel(population)]); 
    mutate_x = population{randnum}; %从种群中随机选择到的一个粒子
    
    
    %对待变异粒子应用变异算子得到新粒子
    x = AddDelete(mutate_x,ub,lb); 
    obj_x = problem.CalObj(x); %计算粒子的目标值
    

    %对种群中的粒子进行筛选,优胜略汰
    cd = 0; 
    cdd = singlenobj;
   
    dominated_index = zeros(size(population,2),1);
    
    for i = 1:numel(population)  %遍历种群中的所有粒子看是否能够支配
        if x(end) == population{i}(end) &&...  %必须保证终点相同才能比较
            boxweakdominates(obj_x(1+singlenobj:singlenobj+singlenobj),obj_population{i}(1+singlenobj:singlenobj+singlenobj),r) %被新粒子在所有方向上盒支配的粒子，会被删除
            dominated_index(i) = 1;
        elseif x(end) == population{i}(end) &&...
            (dominates(obj_population{i}(1+cd:singlenobj+cd), obj_x(1+cd:singlenobj+cd),r) || ...  %至少在一个方向上满足这个条件，就说明他不能加入种群
                dominates(obj_population{i}(1+cdd:singlenobj+cdd), obj_x(1+cdd:singlenobj+cdd),r))
            dominated_index(i,:) = -1;  %若满足则说明新粒子不能添加进种群，有种群中的粒子能够支配或者盒支配他
        end
    end
    % 值为0：终点不同 / 互不支配
    % 值为1：支配掉种群同终点的粒子了
    % 值为-1：被种群中同终点的粒子支配掉了
    
    if ~any(dominated_index == -1)
        population(dominated_index == 1) = []; 
        obj_population(dominated_index == 1) = []; 
        
        population{end+1} = x;  
        obj_population{end+1} = obj_x;  
    end
    

end
epsilon_min = find_min_epsilon(cell2mat(obj_population'), obj_common, nobj, population, common);
final_epsilon = epsilon_min-1;
end


%% 
function mean_IGD = meanIGD(obj_P, obj_common, P, common ,problem)
    min_IGD = zeros(1, size(obj_P, 1));
    for j = 1:size(obj_P, 1)
        matching_indices = find(common(:, end) == P{j}(end));
        matching_common = obj_common(matching_indices, :);   
        min_IGD(j) = MPIGD(obj_P(j,:),problem,matching_common);
    end
    mean_IGD = mean(min_IGD);

end

%% 寻找一定范围内能够满足的最小epsilon
function epsilon_min = find_min_epsilon(obj_P, obj_common, nobj, P, common)

    min_epsilons = zeros(1, size(obj_P, 1));
    
    % 开始二分搜索
    for j = 1:size(obj_P, 1)
        
    % 定义初始范围和精度  近似程度是1+epsilon
        lower_bound = 1+0.01;
        upper_bound = 1+50;
        tolerance = 0.01; % 精度为两位小数
   
        while (upper_bound - lower_bound) > tolerance
            epsilonn = (lower_bound + upper_bound) / 2; % 计算中间值
            result = check(obj_P(j,:), obj_common, epsilonn, nobj, P{j}, common); % 调用目标函数

            if result == 1
                % 如果输出为1，则缩小上界
                upper_bound = epsilonn;
            else
                % 如果输出为0，则提高下界
                lower_bound = epsilonn;
            end
        end
        min_epsilons(j) = round((lower_bound + upper_bound) / 2, 2);
    end
    epsilon_min = mean(min_epsilons);
end


function finalresult = check(obj_P,obj_common,epsilon,nobj,P,common) 
   
    matching_indices = find(common(:, end) == P(end));
    matching_common = obj_common(matching_indices, :);

    result = zeros(1, size(matching_common, 1));


    for i = 1:size(matching_common, 1)

        if epsilondominates(obj_P, matching_common(i, 1:nobj),epsilon)  
            result(i) = true; 
            break;  
        end

    end
    
    finalresult = min(result); 
                              
end


%% 变异操作
function xx = AddDelete(x,ub,lb)  
% Add: 在一个点后面加上一个点
% Delete: 删除一个点后面的点

if rand()<0.5  
    operation = 1;  %添加
else
    operation = 2;  %删除
end

possible_elements = 1:ub; 
if all(ismember(possible_elements, x)) 
    operation = 2;  
end

xx = x; 
mutation_index = randi(length(x));  % 随机选择位置

if operation == 1  %插入
    new_element = randi([lb, ub]); 
    while ismember(new_element, x) 
        new_element = randi([lb, ub]); 
    end
    xx = [x(1:mutation_index), new_element, x(mutation_index+1:end)];
elseif length(x)>2 && operation == 2 && mutation_index < length(x)-1  %删除，删除的时候没有什么约束，直接删除就行了
    xx = [xx(1:mutation_index),x(mutation_index+2:end)];
elseif length(x)>2 && operation == 2 && mutation_index == length(x)-1 %需要删除并且是倒数第二个，那就直接把最后一个位置删除掉
    xx = x(1:mutation_index);
elseif length(x)>2 && operation == 2 && mutation_index == length(x) %需要删除并且是最后一个，那就直接把最后一个位置删除
    xx = x(1:mutation_index-1);
elseif length(x)>2 && operation == 2 && mutation_index == 1 % 需要删除但是是起点，这个时候不进行删除操作
    xx = x;
end

if xx(1) ~= 1 %再去确保路径的起点是不是1
    xx = [1, xx];  % 在x的前面添加1
end


%去重并保持原来顺序
[~, idx] = unique(xx, 'stable');
xx = xx(idx);

end


%% 支配关系的定义
function is_dominated = dominates(x, y, r)
    is_dominated = false;
    rr = r^(length(y));
    if (all(x <= y) && any(x < y) ) ||  (  all(floor(log(x)/log(rr)) <= floor(log(y)/log(rr))) && any(floor(log(x)/log(rr)) < floor(log(y)/log(rr))) )
        is_dominated = true;
    end
end

function is_boxdominated = boxweakdominates(x, y, r) 
    is_boxdominated = false;
    rr = r^(length(y));
    if all(floor(log(x)/log(rr)) <= floor(log(y)/log(rr)))
        is_boxdominated = true;
    end
end

function weak_dominated = epsilondominates(x, y, epsilon)
    weak_dominated = false;
    yy = epsilon * y; %将y的每个元素都进行放大
    if all(x <= yy) 
        weak_dominated = true;
    end
end

