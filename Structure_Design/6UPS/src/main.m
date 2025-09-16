clear;clc;close all
% 本文件用于设计新的6UPS(Stewart)机构
% 机构具体要求见pptx文件
% 开始日期：2025.9.3
% 作者：Lmoyi8

%% 任务拆分
% 1、运动学函数实现：正逆运动学、雅可比矩阵（采用POE方法）
% 2、设计指标分析：
%   （1）静负载
%   （2）行程【优先考虑】
%   （3）刚度
%   （4）分辨率
%   （5）重复定位精度
% 3、优化函数实现：计算评价指标（条件数、LTI）、设计优化函数（初步采用粒子群优化）
% 4、提升点：
%   （1）展示工作空间：方法和形式需要调研
%   （2）尺寸选择

%% 运动学函数
clear;clc;close all
n = 1; % 表示输入输出均为mm
% 采用 Stewart 类中函数
ra =300;%300;
rb = 200;%125;
L = 400;
Robot = Stewart(ra,rb,L,n);

% 静力学检查
% F = rand(3,1)*100;%[100;10;100];
% T = rand(3,1)*100;%[100,0,0]';
% 
% g = Robot.g0;
% 
% qF_force = Robot.force(g,F,T)
% Jac = Robot.Jacobian(g);
% qF_Jac  = Jac'*[F;T]
% 设计指标
% Ra,Rb,H
    x_lim = [-5,5]';
    y_lim = [-5,5]';
    z_lim = [-10,10]'+Robot.H;
    rx_lim = [-10,10]'; % deg
    ry_lim = [-10,10]';
    rz_lim = [-10,10]';
    R = rotz(-10,"deg")*roty(-10,"deg")*rotx(-10,"deg");
    % R = eye(3);
    h0 = Robot.H;
    lo = Robot.L;
    p = [5,-5,h0+10]';
    g = [R,p;0,0,0,1];
    g = Robot.g0;
    s = Robot.ikine(g);
    g = Robot.fkine(s');
    J = Robot.Jacobian(g);
    cond(J)
    Robot.Plot_Robot(g);
    
    % Robot.Robot_Check([x_lim,y_lim,z_lim,rx_lim,ry_lim,rz_lim],1,[500,500,4000,0,0,0]')

%% 优化方法 - 粒子群优化
    
    pop=40; % 种群数量
    dim = 3; % Ra, Rb, H      % 变量维数
    
    
    ub = [300,300,600]; % 变量上界
    lb = [150,100,400]; % 变量下界

    vmax = (ub-lb)/10;%[10,10,10,10]; % 最大速度
    vmin = -vmax; % 最小速度
    maxIter = 40; % 最大迭代次数   

    % fobj = @(X,type)fun(X,type);
    tic;

    [Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,vmax,vmin,maxIter);
    t =toc;

