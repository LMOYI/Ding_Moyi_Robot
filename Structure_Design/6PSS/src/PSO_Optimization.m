%% PSO 粒子群优化
%pop——种群数量
%dim——问题维度
%ub——变量上界，[1,dim]矩阵
%lb——变量下界，[1,dim]矩阵
%fobj——适应度函数（指针）
%MaxIter——最大迭代次数
%Best_Pos——x的最佳值
%Best_Score——最优适应度
clear;clc;close all
%A_L = 400;
pop=100; % 种群数量

dim = 6; % A_d, B_L, B_d, Hy, Da, Db      % 变量维数

% A_d, B_L, B_d, Hy, Da, Db
% a1 - 2点位置角, a2 - 短边圆心角一半
% [a1,a2,b1,b2] in rad

ub = [[ 0 , 30, 0 , 30], 400, 400]; % 变量上界
lb = [[-40, 10,-40, 10], 300, 200]; % 变量下界
vmax = (ub-lb)/10;%[10,10,10,10]; % 最大速度
vmin = -vmax; % 最小速度
maxIter = 40; % 最大迭代次数            
fobj = @(X,type)fun(X,type);
tic;
% -3 同向 -5 反向
type = -3;
[Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,fobj,vmax,vmin,maxIter,type);
t =toc;

%%
a1 = Best_Pos(1);
a2 = Best_Pos(2);
b1 = Best_Pos(3);
b2 = Best_Pos(4);
Da = Best_Pos(5);
Db = Best_Pos(6);
A_L = D;
%% 绘制 PSO 优化结果
figure
plot(IterCurve,'r','linewidth',2);
grid on;
disp(['求解得到的A_d, B_L, B_d, Hy是:',num2str(Best_Pos(1)),' ',num2str(Best_Pos(2)),' ',num2str(Best_Pos(3)),' ',num2str(Best_Pos(4))]);
disp(['最优解对应的函数:',num2str(Best_fitness)]);
 
%% function addition
function [X]=initialization(pop,ub,lb,dim,flag)
    X = zeros(pop,dim);
    for i=1:pop % 逐个种群生成
        for j=1:dim % 逐个变量生成
            X(i,j)=(ub(j)-lb(j))*rand()+lb(j);%在限定的  
        end

        % [a1,a2,b1,b2] in rad

        if flag == 0 % 仅对位置订正，速度不修正
            while abs(X(i,1)-2*X(i,2)) > 90 % 3点位置角超限
                X(i,2)=(ub(2)-lb(2))*rand()+lb(2);% 在限定的  
            end

            while abs(X(i,3)-2*X(i,4)) > 90 % 3点位置角超限
                X(i,4)=(ub(4)-lb(4))*rand()+lb(4);% 在限定的  
            end
        end
    end
end
 
function fitness=fun(X,type)
    % Hy = 150;              % Initial height
    n = 1e3;
    basicX = X(1:4);
    Da = X(5);
    Db = X(6);
    Hy = 120;

    p0 = [0,-Hy,0]'/n;
    R0 = eye(3);
    % fitness = 0;
    if type == -2 || type == -1 || type == -3 || type == -5
        [a,b]=Initial_Config(basicX,type,n,Da,Db);
    else
        [a,b]=Initial_Config(basicX,type,n);
    end
    Len = vecnorm(repmat(p0,1,6) + R0*b - a);
    
    % Len = zeros(1,6);
    % for i = 1:6
    %     Len(i) = norm(R0*b(:,i)+p0-a(:,i)); % 支链长度，单位与n保持一致
    % end

    % 定义参数范围 - 工作空间
    x_values = -25:5:25;
    y_values = 0;
    z_values = -25:5:25;  
    ry_values = -3:2:3;
    rz_values = -3:2:3;
    rx_values = -10:10:10;

    % 生成所有参数组合（使用ndgrid展开）
    [X1, Y, Z, RY, RZ, RX] = ndgrid(x_values, y_values, z_values, ry_values, rz_values, rx_values);

    % 转换为列向量
    params = [X1(:), Y(:), Z(:), RY(:), RZ(:), RX(:)];

    % 初始化存储结果的数组
    num_params = numel(params)/6;  % 总参数组合数
    fitness_values = zeros(num_params, 1);
    s_value = zeros(num_params,6);
    % 启动并行池（如果尚未启动）
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    sdot = @(s1,s2)sum(s1(1:3,:).*s2(4:6,:)+s1(4:6,:).*s2(1:3,:),1); % 定义互易积
    % 并行计算
    parfor i = 1:num_params
        % 解包参数
        x = params(i,1);
        y = params(i,2);
        z = params(i,3);
        ry = params(i,4);
        rz = params(i,5);
        rx = params(i,6);

        % 计算变换矩阵
        pk = [x; y; z]/n;
        Rk = rotz(rz, "deg") * roty(ry, "deg") * rotx(rx, "deg");
        sk = ikine(a,b,Len,pk,Rk); % 计算行程
        s_value(i,:) = sk;
        if ~isreal(sk)
            disp([x,y,z,rx,ry,rz])
            disp(X)
        end
        ak = a;
        ak(2,:) = sk;
        Lk = repmat(pk,1,6) + Rk*b - ak;
        
        % 输入传递指标
        C_phik = ([0,-1,0]*Lk)./vecnorm(Lk,2,1);
        lambda_k = abs(C_phik);

        % 计算雅可比矩阵及其条件数
        [Jac,~] = iJacobian(a, b, Len, pk, Rk); % 计算几何雅可比矩阵，输入单位均一致。
        
        % 输出传递指标
        SO_k = Jac + [skew(pk)*Jac(4:6,:);zeros(3,6)];
        SO_k = SO_k./vecnorm(SO_k(4:6,:),2,1);
        ST_k = [zeros(3,6);Lk]./vecnorm(Lk,2,1);
        ita_k = abs(sdot(SO_k,ST_k))./ (vecnorm(SO_k(1:3,:),2,1).*vecnorm(ST_k(4:6,:),2,1));

        fitness_values(i) = min([ita_k,lambda_k]);
    end

    % 获取全局最优解：LTI = min{Xi}越大越好，算法可求出最小值，故取LTI倒数
    [fitness1, ~] = max(1./fitness_values);
    fitness2 = max(max(s_value) - min(s_value)); % mm
    
    fitness = sqrt(fitness1)*fitness2;
    if isempty(fitness)
        disp('something wrong')
        fitness = 0;
    end
end

 
function [X]=BoundaryCheck(X,ub,lb,dim,flag)
    for i=1:dim
        if X(i)>ub(i)
            X(i)=ub(i);
        end

        if X(i)<lb(i)
            X(i)=lb(i);
        end
    end
    
    if flag == 0 % 仅对位置订正，速度不修正
        while abs(X(1)-2*X(2)) > 90 % 3点位置角超限
            X(2)=(ub(2)-lb(2))*rand()+lb(2);% 在限定的
        end

        while abs(X(3)-2*X(4)) > 90 % 3点位置角超限
            X(4)=(ub(4)-lb(4))*rand()+lb(4);% 在限定的
        end
    end

end
 
function [Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,fobj,vmax,vmin,maxIter,type)
    c1=2.0;
    c2=2.0;
    w_max = 2;
    w_min = 0.4;
    V=initialization(pop,vmax,vmin,dim,1);
    X=initialization(pop,ub,lb,dim,0);

    fitness=zeros(1,pop);

    for i=1:pop
        fitness(i)=fobj(X(i,:),type);
    end
    
    pBest=X;
    pBestFitness=fitness;

    [~,index]=min(fitness);
    gBestFitness=fitness(index);
    gBest=X(index,:);
    
    Xnew=X;
    fitnessNew=fitness;

    for t=1:maxIter
        w = w_max - (w_max - w_min)*t/maxIter;
        
        for i=1:pop
            r1=rand(1,dim);
            r2=rand(1,dim);
            V(i,:)=w*V(i,:)+c1.*r1.*(pBest(i,:)-X(i,:))+c2.*r2.*(gBest-X(i,:));
            V(i,:)=BoundaryCheck(V(i,:),vmax,vmin,dim,1);
            Xnew(i,:)=X(i,:)+V(i,:);
            if type == -2
                if Xnew(6)<Xnew(2)
                    continue
                end
            end
            Xnew(i,:)=BoundaryCheck(Xnew(i,:),ub,lb,dim,0); 

            fitnessNew(i)=fobj(Xnew(i,:),type);
            if fitnessNew(i)<pBestFitness(i)
                pBest(i,:)=Xnew(i,:);
                pBestFitness(i)=fitnessNew(i);
            end

            if fitnessNew(i)<gBestFitness
                gBestFitness=fitnessNew(i);
                gBest=Xnew(i,:);
            end
        end

        X=Xnew;
        fitness=fitnessNew;
        Best_Pos=gBest;
        Best_fitness=gBestFitness;
        IterCurve(t)=gBestFitness;
    end
end