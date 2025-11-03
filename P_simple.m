% 简单算法的截止时间是两方种群都求完自己的PF的时候

% clear
% clc
% close all %清理橱窗
% rng(1) %确保随机数生成是可重复的
% addpath(genpath(pwd)); %将所有目录和子目录添加到运行文件夹下

function t = P_simple(D)
% the hole framework begin
dm = 2; % 参与者数量，一共有两个人
nobj = 4; %目标数量,一共有4个目标
singlenobj = 2; % 每个参与方的有两个目标
common_pf = ones(1,D); %公共解是全1的向量
common_obj = BPAOAZ(common_pf);
common = [common_pf common_obj];

%第一方的pf
first_col = (0:D/2)'; %从0到n/2
second_col = (D:-1:D/2)'; %从n/2到n
pf_1 = [first_col, second_col];

%第二方的pf
first_col = (D:-1:D/2)'; %从 n 到 n/2
second_col = (0:D/2)'; %从 0 到 n/2
pf_2 = [first_col, second_col];

% initialize the population P_1,P_2,Phi
X0 = zeros(1,D); %使用全0向量初始化种群
Y0 = BPAOAZ(X0); % objectives of X0

p_1 = [X0 Y0]; % 初始种群
p_2 = p_1;
phi = p_1;
% begin iterations to uodate population
tic; % 开始计时
while ~( check_1(p_1, D, singlenobj, pf_1) && check_2(p_2, D, singlenobj, pf_2) )  % 简单算法的截止时间是两方种群都求完自己的PF的时候
    
    for p = [1,2] 
        
        %确定本轮遍历的种群以及优化方向
        if p==1
            population = p_1;
            choosed_direction = 0; % 此时优化方向为参与方1
        elseif p==2
            population = p_2;
            choosed_direction = singlenobj; % 此时优化方向为参与方2
        end
        cd = singlenobj - choosed_direction; %用于phi更新的那里
        
        %随机均匀地从种群中选择一个粒子
        randnum = randi([1, size(population,1)]); % 从1到种群粒子数量之间随机选择一个整数
        mutate_x = population(randnum,1:D); %从种群中随机选择到的一个粒子
        x_prime = mutate_x;%待变异的粒子


        %粒子进行进化 one-bit mutation，即随机选择一位然后将其变异
        mutation_index = randi([1, D]); % 随机选择一个位置
        x_prime(mutation_index) = 1-x_prime(mutation_index); %完成一位翻转
        obj_x = BPAOAZ(x_prime); %计算粒子的目标值
        x=[x_prime, obj_x];

        
        
        %对种群中的粒子进行筛选,优胜略汰
        %首先判断这个方向上有没有好的粒子出现，如果没有就不用进行后续判断了
        dominated = false;
        for z = 1:size(population,1) %这里判断的时候用弱支配
            if weakdominates(population(z,D+1+choosed_direction:D+singlenobj+choosed_direction), obj_x(1+choosed_direction:singlenobj+choosed_direction))
                dominated = true;%这时就说明种群中存在比他好或者跟他一样的粒子，这时就不用添加新粒子了。
                break;
            end
        end

        if ~dominated %说明种群中没有比他好的，可以更新种群，去掉种群中被他支配的粒子
            population = [population(~arrayfun(@(z) dominates(... %有相等的粒子的话在上一步就被卡掉了。能到这一步说明没有相等的，用大于等于没事
                obj_x(1+choosed_direction:singlenobj+choosed_direction), population(z,D+1+choosed_direction:D+singlenobj+choosed_direction)), 1:size(population,1)), :); x];
            
            
            %更新原本的种群
            if p==1
                p_1 = population;
            elseif p==2
                p_2 = population;
            end
            
            
            %判断phi是否需要更新
            dominated_in_phi = false; 
            for z = 1:size(phi,1)
                if weakdominates(phi(z,D+1:D+nobj), obj_x(1:nobj))
                    dominated_in_phi = true; %这时说明种群中存在粒子能够在所有参与方上弱支配新生成的粒子，这时就不会再对phi进行更新了
                    break;
                end
            end

            if ~dominated_in_phi %更新phi，去掉phi中至少在一个方向上比他差的
                phi = [phi(~arrayfun(@(z) dominates(...  %只要在一方上被支配就删掉他
                    obj_x(1+choosed_direction:singlenobj+choosed_direction), phi(z,D+1+choosed_direction:D+singlenobj+choosed_direction)) || ...
                    dominates(obj_x(1+cd:singlenobj+cd), phi(z,D+1+cd:D+singlenobj+cd)), 1:size(phi,1)), :); x];
            end
        end
     
        
    end
    
end
t = toc; % 结束计时并返回运行时间
% fprintf('simple算法dim=%.0f时运行时间为 %.7f 秒\n', D,t);



%% 检查种群是否满足要求
function isEqual = check_1(P, D, singlenobj, pf) %pf_1
    sortedP = sortrows(P, D+1);
    isEqual = isequal(sortedP(:,D+1:D+singlenobj), pf);
end

function isEqual = check_2(P, D, singlenobj, pf) %pf_2
    sortedP = sortrows(P, D+singlenobj+singlenobj);
    isEqual = isequal(sortedP(:,D+1+singlenobj:D+singlenobj+singlenobj), pf);
end

%% 支配关系的定义
function is_dominated = dominates(x, y) %支配，不包含相等的情况
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

end

