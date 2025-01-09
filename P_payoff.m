% clear
% clc
% close all %清理橱窗
% rng(1) %确保随机数生成是可重复的
% addpath(genpath(pwd)); %将所有目录和子目录添加到运行文件夹下
function t = P_payoff(D)
% the hole framework begin
dm = 2; % 参与者数量
nobj = 4; %目标数量
singlenobj = 2; 
common_pf = ones(1,D); 
common_obj = BPAOAZ(common_pf);
common = [common_pf common_obj];

% initialize the population P
X0 = zeros(1,D); %使用全0向量初始化种群
Y0 = BPAOAZ(X0); % objectives of X0

population = [X0 Y0]; % 初始化种群
% begin iterations to uodate population
tic; % 开始计时
while ~isequal(population, common)
    %随机均匀地从种群P中选择一个粒子
    randnum = randi([1, size(population,1)]);
    mutate_x = population(randnum,1:D); 
    mutate_y = mutate_x; 
    
    %被选择到的粒子进行进化 one-bit mutation，即随机选择一位然后将其变异
    mutation_index = randi([1, D]); 
    mutate_y(mutation_index) = 1 - mutate_y(mutation_index); 
    
    obj_mutate_x = BPAOAZ(mutate_x);
    x =[mutate_x  obj_mutate_x];

    obj_mutate_y = BPAOAZ(mutate_y);%变异之后的
    y =[mutate_y  obj_mutate_y];
    
    %根据总收益决定是否添加新粒子
    MP_payoff = 0; 
    for d = 0:dm-1 
        MP_payoff = MP_payoff + Payoff(x,y,d,D,singlenobj);
    end
    
    %对变异得到的粒子进行筛选
    if MP_payoff >= 0 
        population(randnum, :) = []; 
        population = [population; y];
    end
 
end

t = toc; % 结束计时并返回运行时间




%% 收益的定义
function payoff = Payoff(x,y,d,dim,singlenobj) % 求的是从x进化到y的收益
    x_obj = x(dim+1+singlenobj*d : dim+singlenobj+singlenobj*d);
    y_obj = y(dim+1+singlenobj*d : dim+singlenobj+singlenobj*d);
    if all(x_obj <= y_obj) 
        payoff = sum(y_obj - x_obj); 
    elseif all(x_obj >= y_obj) 
        payoff = sum(y_obj - x_obj); 
    else % 互不支配
        payoff = 0;
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
