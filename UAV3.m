% 总顶点数量为30
classdef UAV3 < handle
        %% Properties  
        properties (SetAccess = immutable)
            M;         
            DM;
            data;
            lower;
            upper;  
            D; % 让他固定维度，否则没办法将新的粒子直接添加进种群
            encoding;
            funname;
            %maxFE;
            seed;
            points;
            epsilon;
            lb;
            ub;
            common;
            wheight;
        end
        
        properties (SetAccess = private) %私有属性，用来记录评估次数
            calcount;
        end     

        methods
        %% Initialization
        function obj = UAV3(seed)
            obj.seed = seed; %设置随机数种子，确保每次实验的可重复性。
            assert(obj.seed<=30 && obj.seed>=1, "The random seed must be an integer in the interval [1,30]");
            
            %obj.maxFE = 10000;
            obj.calcount = 0; %计算函数的评估次数
            
            obj.points = 30;  %总顶点数为30
            
            obj.M = 4; %总的目标数量
            obj.DM = 2; %参与方数量
            obj.data = load('data.mat');
            obj.data = obj.data.data; %加载文件中的具体数据
            
            % 这里的上界和下界是用来控制轨迹点的范围
            obj.lb = 2; %为避免去重和去环操作过于复杂，这里使最小的顶点直接为2
            obj.epsilon = 1+0.1; % 近似程度
            obj.ub = obj.points; %最大的顶点序号
            
            obj.D = obj.ub - 1; %这个是简单路径的长度，所以就要去重、去环，才能达到这个长度
            obj.encoding = 'real'; %解是通过实数（连续值）表示的。还是连续值，只是后面给他取整得到轨迹路径了
            obj.funname = {'length','distance','fatal','eco'}; %四个目标函数的名称
            % 这里的上界和下界是用来控制高度的
            obj.lower = 0; %上界和下界，用于限制高度的范围，值是根据data中的lb来的
            obj.upper = 1;  %上界和下界，用于限制高度的范围，值是根据data中的ub来的
            
            obj.wheight = 0.1; %直接设置权重，用于求高度的目标值那里
            obj.common = obj.data.common3; %顶点数量为30时的公共解
        end

        %% Calculate population's objective values for each party
        function PopObj = CalObj(obj,PopDec)         
            PopObj = zeros(size(PopDec,1),obj.M); %初始化目标值矩阵，粒子个数*目标函数数量
%             assert(obj.calcount<=obj.maxFE, "Maximum number of evaluations exceeded.");%如果评估次数不够就不计算了

            for i =1:size(PopDec,1)
                PopObj(i,:)=obj.Calobj(PopDec(i,:),obj.data);%这个是计算单个个体的目标值
            end
            obj.calcount = obj.calcount + size(PopDec,1);
        end

        %% Calculate single individual's objective values for each party     
        function [fit,trace_ijz]=Calobj(obj,x,data) 
            data = obj.data;
            trace_code = x ;%变量x中每个元素进行“取整”操作。再根据前面生成粒子的时候，所以后面就不需要在进行大小范围的限制了
            trace=[[1:length(trace_code)]' trace_code'];
            
            %长度越小权重越小，这是再让他往公共解那里靠
            height_code = length(trace_code)-1;   %在这里设置的就是一条路径一个值
            %height_code = min(max(height_code, obj.lower), obj.upper); %进行范围的限制
                    
            %初始化各种变量，开始计算目标值
            trace_heightmin=zeros(1,length(trace));
            trace_height=zeros(1,length(trace));
            risk_die = zeros(1,length(trace));
            risk_money = zeros(1,length(trace));
            
            %disp(trace_code);
            for i =1:length(trace) %遍历长度/轨迹点的每一个元素
                %trace_heightmin(i)=max(data.minh(trace(i,1),trace(i,2)),5); %每个轨迹点的最小高度
                %trace_height(i)=trace_heightmin(i)+(data.maxh-trace_heightmin(i))*height_code; %每个轨迹点的实际高度
                h = height_code;
                
                % 根据UAV高度和人群分布populations_risk计算的死亡风险
                %Cpf=getC_Risk(getR_pf(getV(h,data),data),data.populations_risk(trace(i,1),trace(i,2)),data); %计算风险
                %根据路径点的道路风险road_risk计算的死亡风险
                Cvf=getC_Risk(data.R_vf,data.road_risk(trace(i,1),trace(i,2)),data); %计算风险
                %risk_die(i)=Cpf+Cvf; %每个路径点的总死亡风险
                risk_die(i)=Cvf;
                
                C_rpd=getC_rpd(h,data); %计算路径高度h下的UAV财务成本，比如飞行设备损耗和维护费用
                risk_money(i)=C_rpd; %每个路径点的财务风险
            end
            trace_xyz=[(trace-1).*data.map_step trace_height']; %计算轨迹在三维空间下的坐标，然后根据三维空间下的坐标去计算目标值
            
            %计算轨迹之间的夹角，表示路径的变化
            
            %根据坐标求夹角，路径长度啥的
            trace_ijz=[trace trace_height'];
            trace_len = sqrt(sum((trace_xyz(1:end-1,:)-trace_xyz(2:end,:)).^2,2)); % 路径中每一段线段的长度，计算的是三位欧几里得距离
            
            %目标函数值计算
            disMat = min(pdist2(data.IOT_pos*100,trace_xyz),[],2); 
            %计算路径点和设备位置(需要监控或覆盖的区域)之间的距离。对于每个设备，选择与路径点之间的最短距离
            distsum =  sum(disMat); %距离——目标2，对所有设备的最短距离求和
            
            %distsum =  length(trace_code); %单源最短路径，经过的点越少越好，这样能在路上浪费更少的时间
            tracelensum = sum(trace_len);  %轨迹长度——目标1，多有路径段长度的总和，表示UAV飞行路径的总距离
            
            risk_diesum = sum(risk_die);%死亡风险——目标3，所有路径点死亡风险的总和
            risk_moneysum = sum(risk_money); %财务风险——目标4，所有路径点财务风险的总和 
            
            %计算目标值：轨迹长度、距离、死亡风险、财务风险
            fit=[tracelensum;distsum;risk_diesum;risk_moneysum];
         end
    end
end