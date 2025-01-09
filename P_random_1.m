
function t = P_random_1(D,PHI)
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

population = [X0 Y0]; % 初始种群
% begin iterations to uodate population
tic; % 开始计时
while ~isequal(population, common)
    %随机均匀地从种群P中选择一个粒子
    randnum = randi([1, size(population,1)]); 
    mutate_x = population(randnum,1:D); 
    x = mutate_x;
    
    %该粒子随机选择一个方向进化
    if rand()<PHI %\phi^+ 暂时设置为0.5
        choosed_direction = 0; 
    else
        choosed_direction = singlenobj; 
    end
    
    %选择完方向之后，粒子进行进化 one-bit mutation，即随机选择一位然后将其变异
    mutation_index = randi([1, D]); 
    x(mutation_index) = 1-x(mutation_index); 
    obj_x = BPAOAZ(x); 

    
    %对种群中的粒子进行筛选,优胜略汰
    dominated_index = zeros(size(population, 1), 1); 
    for i = 1:size(population, 1)  
        if dominates(obj_x(1+choosed_direction:singlenobj+choosed_direction),population(i,D+1+choosed_direction:D+singlenobj+choosed_direction))
            dominated_index(i,:) = 1; 
        elseif weakdominates(population(i,D+1+choosed_direction:D+singlenobj+choosed_direction),obj_x(1+choosed_direction:singlenobj+choosed_direction))
            dominated_index(i,:) = -1; 
        else
            dominated_index(i,:) = 0; 
        end
    end
    
    if ~any(dominated_index == -1)
        population(dominated_index == 1, :) = []; 
        xs =[x  obj_x];
        population = [population;xs];
    end
    
    
    %最后只保留种群中在这个方向上的非支配解
    if size(population,1) > 1
        population = removeNondominated(population,choosed_direction,singlenobj,D);
    end

end
t = toc;


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

function newpopulation = removeNondominated(population,d,singlenobj,D) %删除种群中在这个方向上的非支配解
    n = size(population, 1); % 种群大小
    is_dominated = false(n, 1); % 标记每个解是否被支配

    % 遍历种群中的每一对粒子
    for i = 1:n
        for j = 1:n %对于第i个粒子，判断他是不是被种群中的其他粒子支配
            if i ~= j 
                % 如果粒子j支配粒子i，则标记i为被支配的粒子
                if weakdominates(population(j,D+1+d:D+singlenobj+d), population(i,D+1+d:D+singlenobj+d))
                    is_dominated(i) = true;
                    break; % 一旦被支配，无需再检查
                end
            end
        end
    end

    % 过滤出非支配解
    newpopulation = population(~is_dominated, :);

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