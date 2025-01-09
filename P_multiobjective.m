

% clear
% clc
% close all %清理橱窗
% rng(1) %确保随机数生成是可重复的
% addpath(genpath(pwd)); %将所有目录和子目录添加到运行文件夹下

function t = P_multiobjective(D)
% the hole framework begin
nobj = 4; %目标数量,一共有4个目标

% initialize the population P
X0 = zeros(1,D); %使用全0向量初始化种群
Y0 = AOAZ(X0); % objectives of X0

population = [X0 Y0]; % 初始种群

%四目标的pf
pf = generatepf(D,nobj);


tic; % 开始计时
while ~check(population, D, nobj, pf) 

    %随机均匀地从种群中选择一个粒子
    randnum = randi([1, size(population,1)]); 
    mutate_x = population(randnum,1:D); 
    x_prime = mutate_x;


    %粒子进行进化 one-bit mutation，即随机选择一位然后将其变异
    mutation_index = randi([1, D]); 
    x_prime(mutation_index) = 1-x_prime(mutation_index); 
    obj_x = AOAZ(x_prime); 
    x=[x_prime, obj_x];



    %对种群中的粒子进行筛选,优胜略汰
    dominated = false;
    for z = 1:size(population,1)
        if weakdominates(population(z,D+1:D+nobj), obj_x(1:nobj))
            dominated = true;
            break;
        end
    end

    if ~dominated 
        population = [population(~arrayfun(@(z) dominates(... 
            obj_x(1:nobj), population(z,D+1:D+nobj)), 1:size(population,1)), :); x];
    end
   
end
t = toc; 



%% 生成部分PS，得到四目标的PF
function pf = generatepf(D,nobj)
    pf = zeros(D/2+1,nobj);

    p = ones(D/2, D/2); %一半
    pp = ones(D/2, D/2); %一半
    x = ones(1,D);
    for i = D/2:-1:1  %1的数量不断增多
        pp(D/2+1-i, 1:i) = 0;  % 前 i-1 个元素设为 0
    end
    
    p1 = [p, pp]; %第一方的
    p2 = [pp, p]; %第二方的
    
    ppp = [p1;p2;x];
    
    for j = 1:size(ppp,1)
        pf(j,1:nobj) = AOAZ(ppp(j,:));  
    end
end


%% 检查种群是否满足要求
function isEqual = check(P, D, nobj, pf) %pf
    sortedP = sortrows(P, [D+1, D+nobj]);
    isEqual = isequal(sortedP(:,D+1:D+nobj), pf);
end


%% 支配关系的定义,该问题为最大化问题
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
function result = AOAZ(x)
    % AOAZ: Computes the pseudo-Boolean function BPAOAZ.
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

