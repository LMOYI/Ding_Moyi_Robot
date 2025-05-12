%6PSS_COMPUTE 建立6PSS运动学模型，并计算指定负载下静力学
% 参考文献：
%   [1] 沈音诚学长硕士论文
%   [2] X.L Yang, H.T Wu, etc. A dual quaternion solution to the forward kinematics of a class of six-DOF parallel robots with full or reductant actuation,
%   [3] 一篇法文的关于雅可比矩阵的文章
% 调用函数：
% 机构基本信息：
%   1、Dof = 6；采用ZYX欧拉角；
%   2、负载：90kg
%   3、XY平动范围：+-25mm；三方向转动范围：+-3°
%   4、推杆速度：10mm/s；推杆加速度：10mm/s^2；连杆长度：l>100mm;
clear;clc;close all;
n =1e3;              % control the dimension: n=1 is in mm, n=1000 is in m
Hy = 150;              % Initial height
p0 = [0,-Hy,0]'/n;
R0 = eye(3);
%下述值需要被优化，含义见ppt
A_L = 400;          % Distance btw lines of the Base
A_d = 80;%437.7116;          % Distance btw pts on each line of the Base

B_L = 400;%453.5131;%B_d*ratio;          % Distance btw lines of the Platform
B_d = 200;%616.6242;          % Distance bte pts on each line of the Platform


%初始平台布局生成，调用函数 [a,b]=Initial_Config(A_d, A_L, B_d, B_L,k,n)
[a,b]=Initial_Config(A_d, A_L, B_d, B_L,-1,n,500,500);
Len = vecnorm(repmat(p0,1,6)+R0*b - a);
% Len = zeros(1,6);
% for i = 1:6
%     Len(i) = norm(R0*b(:,i)+p0-a(:,i));
% end
% robot_plot(a, b, zeros(6,1), p0, R0)
% plot(a(1,:),a(3,:),'r-')
% hold on
% plot(b(1,:),b(3,:),'b--')
% legend('基座','动平台')


%
% clc
skew=@(x)[0,-x(3),x(2); x(3),0,-x(1);-x(2),x(1),0];
R1 = rotz(-30, "deg") * roty(3, "deg") * rotx(3, "deg");
p1 = [-25,0,-25]'/n;
s0 = ikine(a,b,Len,p0,R0);
a0 = a;
a0(2,:) = s0;
L_check = vecnorm(repmat(p0,1,6) + R0*b - a0 ,2,1);

robot_plot(a,b,s0,p0,R0)

% sk = [10,10,10,11,11,11]/n;
pk = p0+[0,0,0]'/n;
x = 0;      
angle_x = x:0.1:x+pi/12;
pic_num = length(angle_x);
cond_max = 0;
k = 0;
    % 定义参数范围
    x_values = -25:5:25;
    % y_values = -25:5:25;
    y_values = 0;
    z_values = -25:5:25;  
    ry_values = -3:2:3;
    rx_values = -3:2:3;
    rz_values = -30:10:30;

    % 生成所有参数组合（使用ndgrid展开）
    [X, Y, Z, RY, RZ, RX] = ndgrid(x_values, y_values, z_values, ry_values, rz_values, rx_values);

    % 转换为列向量
    params = [X(:), Y(:), Z(:), RY(:), RZ(:), RX(:)];

    % 初始化存储结果的数组
    num_params = numel(params)/6;  % 总参数组合数
    fitness_values = zeros(num_params, 1);
    fitness3 = zeros(num_params, 1);

    % 启动并行池（如果尚未启动）
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    sdot = @(s1,s2)sum(s1(1:3,:).*s2(4:6,:)+s1(4:6,:).*s2(1:3,:),1);
    % 并行计算
    parfor i = 1:num_params
    % i=1;
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

        % 计算雅可比矩阵及其条件数
        sk = ikine(a,b,Len,pk,Rk)
        ak = a;
        ak(2,:) = sk;
        Lk = repmat(pk,1,6) + Rk*b - ak; % 传力方向
        C_phik = ([0,-1,0]*Lk)./vecnorm(Lk,2,1);
        lambda_k = abs(C_phik);


        [Jac,~] = iJacobian(a, b, Len, pk, Rk);
        % SO_k = Jac./vecnorm(Jac(4:6,:),2,1);    % -> matrix<6,num> -> [x,y,z,w1,w2,w3]'
        SO_k = Jac + [skew(pk)*Jac(4:6,:);zeros(3,6)];
        SO_k = SO_k./vecnorm(SO_k(4:6,:),2,1);
        ST_k = [zeros(3,6);Lk];
        % ST_k = Jac_inv'./vecnorm(Jac_inv(:,4:6)',2,1);
        ita_k = abs(sdot(SO_k,ST_k))./ (vecnorm(SO_k(1:3,:),2,1).*vecnorm(ST_k(4:6,:),2,1));

        fitness3(i) = min([ita_k,lambda_k]);
        fitness_values(i) = cond(Jac);
    end

   % 获取全局最优解
    [cond_max, max_i] = max(fitness_values);
    [LTI_min,min_j] = min(fitness3);
    disp('传递系数')
    disp(LTI_min)
    disp('雅可比条件数')
    disp(cond_max)
%%
for i=-pi/6:0.0001:pi/6
    
    Rk = rotz(i);
    sk = ikine(a,b,Len,pk,Rk);
    Jac = iJacobian(a,b,Len,pk,Rk);
    if cond(Jac)>cond_max
        cond_max = cond(Jac);
        k = i;
    end
end

    % if cond(Jac)>1e5
        % robot_plot(a,b,sk,pk,Rk)
        % i
    % end
    % close all
    % currFrame(pic_count) = getframe;

%% PSO
%pop——种群数量
%dim——问题维度
%ub——变量上界，[1,dim]矩阵
%lb——变量下界，[1,dim]矩阵
%fobj——适应度函数（指针）
%MaxIter——最大迭代次数
%Best_Pos——x的最佳值
%Best_Score——最优适应度
clear;clc;close all
A_L = 400;
pop=50;
dim=4; % A_d, B_L, B_d, Hy 
ub = [300, 400, 300, 200];
lb = [100, 200, 100, 100];
vmax=[25,50,25,10];
vmin=-vmax;
maxIter=20;
fobj=@(X)fun(X);
tic;
[Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,fobj,vmax,vmin,maxIter);
t =toc;
%%
figure
plot(IterCurve,'r','linewidth',2);
grid on;
disp(['求解得到的A_d, B_L, B_d, Hy是:',num2str(Best_Pos(1)),' ',num2str(Best_Pos(2)),' ',num2str(Best_Pos(3)),' ',num2str(Best_Pos(4))]);
disp(['最优解对应的函数:',num2str(Best_fitness)]);
 

%%
R = @(roll,pitch,yaw)rotz(yaw,"deg")*roty(pitch,"deg")*rotx(roll,"deg");
F_joint = zeros(3,6,51*51);
% n1 =[];N=[];
F = [0;-1200;0];
T = [0;0;0];
y = 0; num=0;
num_F = 51*51*7*7*13;
S_length = zeros(6,num_F);
    % 定义参数范围
    % x_values = -25:5:25;
    % y_values = -25:5:25;
    % z_values = -25:5:25;  % 注意：这个范围实际只会生成25
    % ry_values = -3:2:3;
    % rx_values = -3:2:3;
    % rz_values = -30:10:30;

    % % 生成所有参数组合（使用ndgrid展开）
    % [X, Y, Z, RY, RZ, RX] = ndgrid(x_values, y_values, z_values, ry_values, rz_values, rx_values);
    % 
    % % 转换为列向量
    % params = [X(:), Y(:), Z(:), RY(:), RZ(:), RX(:)];

    % 初始化存储结果的数组
    % num_params = numel(params)/6;  % 总参数组合数
    % fitness_values = zeros(num_params, 1);

    % 启动并行池（如果尚未启动）
    % if isempty(gcp('nocreate'))
    %     parpool('local');
    % end

    % 并行计算
    % parfor i = 1:num_params
    %     % 解包参数
    %     x = params(i,1);
    %     y = params(i,2);
    %     z = params(i,3);
    %     ry = params(i,4);
    %     rz = params(i,5);
    %     rx = params(i,6);
    % 
    %     % 计算变换矩阵
    %     pk = [x; y; z];
    %     Rk = rotz(rz, "deg") * roty(ry, "deg") * rotx(rx, "deg");
    % 
    %     % 计算雅可比矩阵及其条件数
    %     Jac = iJacobian(a, b, Len, pk, Rk);
    %     fitness_values(i) = cond(Jac);
    % end

    % 获取全局最优解
    % [CondJ, ~] = max(fitness_values);
    count = 0;
for x= -25:1:25
    for z = -25:1:25
        % N1 = 0;
        Fx = zeros(1,6); Fy = zeros(1,6);Fz = zeros(1,6);
        p = [x;y;z]/n;
        for roll = -3:1:3
            for pitch =-3:1:3
                for yaw = -30:5:30
                    R1=R(roll,pitch,yaw);
                    s = ikine(a,b,Len,p,R1);
                    count = count +1;
                    S_length(:,count) = s;
                    % Jac = iJacobian(a,b,Len,p,R1);
                    % n0 = cond(Jac);
                    % n1 = [n1,n0];
                    % if n0>N1
                    %     N1=n0;
                    % end
                    [f,~] = force(a,b,s,p,R1,F,T);
                    for k = 1:6
                        if abs(f(1,k))>Fx(k)
                            Fx(k) = f(1,k);
                        end
                        if abs(f(2,k))>Fy(k)
                            Fy(k) = f(2,k);
                        end
                        if abs(f(3,k))>Fz(k)
                            Fz(k) = f(3,k);
                        end
                    end
                end
            end
        end
        % N = [N,N1];
        num = num+1;
        F_joint(:,:,num) = [Fx;Fy;Fz];
    end
end
% 
%%
N_plot= reshape(N,51,51);
contour(-25:1:25,-25:1:25,N_plot'');
xlabel("x/mm");
ylabel("y/mm");
title("最大雅可比条件数图");
%%
clc
F1 = reshape(F_joint(:,1,:),3,2601);
F2 = reshape(F_joint(:,2,:),3,2601);
F3 = reshape(F_joint(:,3,:),3,2601);
F4 = reshape(F_joint(:,4,:),3,2601);
F5 = reshape(F_joint(:,5,:),3,2601);
F6 = reshape(F_joint(:,6,:),3,2601);

disp('关节1最大受力')
disp(max(vecnorm(F1)))
disp('关节2最大受力')
disp(max(vecnorm(F2)))
disp('关节3最大受力')
disp(max(vecnorm(F3)))
disp('关节4最大受力')
disp(max(vecnorm(F4)))
disp('关节5最大受力')
disp(max(vecnorm(F5)))
disp('关节6最大受力')
disp(max(vecnorm(F6)))

%% 行程
clc
S = (max(S_length,[],2) - min(S_length,[],2))*n
S_1 = max(S_length,[],2)*n
S_2 = min(S_length,[],2)*n
%%
disp('关节1行程范围')
disp(max(vecnorm(F1)))
disp('关节2最大受力')
disp(max(vecnorm(F2)))
disp('关节3最大受力')
disp(max(vecnorm(F3)))
disp('关节4最大受力')
disp(max(vecnorm(F4)))
disp('关节5最大受力')
disp(max(vecnorm(F5)))
disp('关节6最大受力')
disp(max(vecnorm(F6)))
%%
clc
s = [10,11,10,10,10,10]/n;
[p,R1] = fkine(a,b,Len,s,n);
robot_plot(a,b,s,p/n,R1)

F = [0 0 -1200]';   % 末端载荷
T = [0 0 0]';   % 末端扭矩
[f1,f2] = force(a,b,s,p(1:3)/n,R1,F,T);
Jac = iJacobian(a,b,Len,p/n,R1);
drive_forceref = -(Jac'*[F;T])';
% check

% a1 = R1*b;
% a2 = a; a2(2,:) = s;
% T1 = zeros(3,1);
% for i =1:6
% T1 = T1+skew(a1(:,i))*f1(:,i);
% L = p(1:3)/n + a1(:,i) - a2(:,i);
% skew(L)*f1(:,i)
% end


f1;
disp('雅可比矩阵计算力为：')
disp(drive_forceref)
disp('受力分析计算力为：')
disp(f1)

check = (f1(2,:) - drive_forceref)./f1(2,:);
disp('相对误差为:')
disp(check)
% % a_drive = [10,10,10,10,10,10]/n;
% % [p,R]=fkine(a,b,Len,a_drive,n);
%% 验证雅可比矩阵
clear;clc;close all;

n = 1000;              % control the dimension: n=1 is in mm, n=1000 is in m
Hy = 400;              % Initial height
p0 = [0,-Hy,0]'/n;
R0 = eye(3);
%下述值需要被优，含义见ppt
A_d = 140;          % Distance btw pts on each line of the Base
A_L = 600;          % Distance btw lines of the Base
B_d = 140;          % Distance bte pts on each line of the Platform
B_L = 300;          % Distance btw lines of the Platform

%初始平台布局生成，调用函数 [a,b]=Initial_Config(A_d, A_L, B_d, B_L,k,n)
[a,b]=Initial_Config(A_d, A_L, B_d, B_L,2,n);

Len = zeros(1,6);
for i = 1:6
    Len(i) = norm(R0*b(:,i)+p0-a(:,i));
end
skew=@(x)[0,-x(3),x(2); x(3),0,-x(1);-x(2),x(1),0];
q=@(t)-[30*t;20*t;10*t;-10*t;20*t;0]/n;
dq = -[30;20;10;-10;20;0]/n;
% q = @(t)[t,t,t,t,t,t]';
% dq = [1,1,1,1,1,1]'/n;
dt=1e-3;
t = 0:dt:1;
num = length(t);
v1 = zeros(6,num);
% [p1,R1] = fkine(a,b,Len,q(2),n);
% s = ikine(a,b,Len,p1/n,R1);
% [p2,R2] = fkine(a,b,Len,s,n);

% g2 = zeros(4,4,num);
p2 = zeros(3,num);
R2 = zeros(3,3,num);
for i =1:num
    [p1,R1] = fkine(a,b,Len,q(t(i)),n);
    s1 = ikine(a,b,Len,p1/n,R1);
    Jac = iJacobian(a,b,Len,p1/n,R1);
    v1(:,i) = Jac*dq;
    v1(1:3,i)=v1(1:3,i)*n;
    p2(:,i) = p1(1:3);
    R2(:,:,i) = R1;
end
% g20 =
V2 = zeros(6,num-1);
for i =1:num-1
    V2(1:3,i) = (p2(:,i+1)-p2(:,i))/dt;
    % R_t = R2(:,:,i+1)/R2(:,:,i);
    % xi_hat = logm(R_t);
    xi1 = logm(R2(:,:,i));
    xi2 = logm(R2(:,:,i+1));
    V2(4:6,i) = ([xi2(3,2),xi2(1,3),xi2(2,1)]'-[xi1(3,2),xi1(1,3),xi1(2,1)]')/dt;
    % V2(4:6,i) = [xi_hat(3,2),xi_hat(1,3),xi_hat(2,1)]';
    % g = g2(:,:,i+1)/g2(:,:,i)
    % xi_hat = logm(g);
    % omega = [xi_hat(3,2),xi_hat(1,3),xi_hat(2,1)]';
    % V2(4:6,i) = omega/dt;
    % V2(1:3,i) = (g2(1:3,4,i+1) - g2(1:3,4,i))/dt;
end
v1(:,1:10) % in mm
V2(:,1:10) %

%% forward kinematics
% 对偶四元数方法

clc
% [w,theta] = Quaternion(R);
% sigma = [w.*sind(theta/2);cosd(theta/2)];

% unit is meter
sigma_0  = [0 0 0 1]';      %
lambda_0 = [0 -Hy 0 0]'/n;

% a = [184.6 45.2 0;-53.2 182.4 0;-131.4 137.2 0;-131.4 -137.2 0;-53.2 -182.4 0;184.6 -45.2 0]';
% b = [131.4 137.2 0;53.2 182.4 0;-184.6 45.2 0; -184.6 -45.2 0;53.2 -182.4 0; 131.4 -137.2 0]';
% L = [220,250,250,180,180,190];

Q = zeros(8,8,8);
Q(1:4,1:4,7)=2*eye(4);
Q(5:8,1:4,8)=eye(4);
Q(1:4,5:8,8)=eye(4);

C = [Len.^2,1,0]';

a_drive = [100,200,0,50,50,0];
a(2,:) = a_drive;

for i = 1:6
    Ai = a(:,i)+b(:,i);
    Bi = b(:,i)-a(:,i);
    Qi11 = [2*(Bi(1)^2+Ai(2)^2+Ai(3)^2)     , -4*(a(1,i)*b(2,i)+a(2,i)*b(1,i) ) ,  -4*(a(1,i)*b(3,i)+a(3,i)*b(1,i) ),  -4*(-a(2,i)*b(3,i)+a(3,i)*b(2,i));
        -4*(a(1,i)*b(2,i)+a(2,i)*b(1,i) ),  2*(Ai(1)^2+Bi(2)^2+Ai(3)^2 )     ,  -4*(a(2,i)*b(3,i)+a(3,i)*b(2,i) ),  -4*(-a(1,i)*b(3,i)+a(3,i)*b(1,i) );
        -4*(a(1,i)*b(3,i)+a(3,i)*b(1,i) ), -4*(a(2,i)*b(3,i)+a(3,i)*b(2,i) ) ,   2*(Ai(1)^2+Ai(2)^2+Bi(3)^2 )    ,  -4*(-a(1,i)*b(2,i)+a(2,i)*b(1,i) );
        -4*(-a(2,i)*b(3,i)+a(3,i)*b(2,i)), -4*(-a(1,i)*b(3,i)+a(3,i)*b(1,i) ), -4*(-a(1,i)*b(2,i)+a(2,i)*b(1,i) ),   2*(Bi(1)^2+Bi(2)^2+Bi(3)^2 )     ];
    Qi12 = [0      , 2*Ai(3),-2*Ai(2), 2*Bi(1);
        -2*Ai(3), 0      , 2*Ai(1), 2*Bi(2);
        2*Ai(2),-2*Ai(1), 0      , 2*Bi(3);
        -2*Bi(1),-2*Bi(2),-2*Bi(3), 0]';
    Qi21 = Qi12';
    Qi22 = 2*eye(4);
    Q(:,:,i) = [Qi11,Qi12;Qi21,Qi22];
end

flag = 1;
X = [sigma_0;lambda_0];
% X =
% for i = 1:8
%     flag1 = Q(:,:,i)*X
%     flag2 = Q(:,:,i)-Q(:,:,i)'
% end


% for i = 1:6
% 1/2*X'*Q(:,:,i)'*X-C(i)
% end
% X = [0 0 0 1 0 1 0 0]';

while norm(flag)>1e-10
    J = [Q(:,:,1)*X, Q(:,:,2)*X, Q(:,:,3)*X, Q(:,:,4)*X, Q(:,:,5)*X, Q(:,:,6)*X, Q(:,:,7)*X, Q(:,:,8)*X]';
    % DET_J = det(J)
    % % delta_X = (J'*J)\(J'*C);
    delta_X = pinv(J)*C;
    flag = 0.5*X - delta_X;
    X = 0.5*X + delta_X;
end
sigma  = X(1:4);
lambda = X(5:8);
theta = 2*acosd(sigma(4));
w = sigma(1:3)/sind(theta/2);
p = quaternion_multi(lambda,conjugate(sigma));
R1 = rotation_matrix(w,theta);


%% Inverse kinematics
%Given configuration of end-effector [x,y,z,alpha,beta,gamma]
%Find relative configurations of actuators
alpha = 0;
beta  = 0;
gamma = 0;
p = p0+[0,200,0]'/n;

% 绕定轴ksi 旋转theta
ksi = [0,1,0]';
ksi_hat = skew(ksi);
theta = 0;
% R = Rev(gamma,3)*Rev(beta,2)*Rev(alpha,1);
R1 = expm(ksi_hat*theta);

S = zeros(1,6); % P副上球铰相对位置（在y向）-反映驱动距离

for i = 1:6
    B = p+R1*b(:,i);
    S(i) = B(2)+sqrt(Len(i)^2 - (B(1)-a(1,i))^2 - (B(3)-a(3,i))^2);
end


%% Jacobian
% S^dot = J^inv * V
clc
p = p(1:3);
Jac_inv = zeros(6,6);
for i = 1:6
    g = R1*b(:,i);
    rho = sqrt(Len(i)^2-(p(1)+g(1)-a(1,i))^2-(p(3)+g(3)-a(3,i))^2 ); % 连杆向量y向分量
    r =  [p(1)+g(1)-a(1,i), -rho, p(3)+g(3)-a(3,i)]';   %
    G = -skew(g);
    Jac_inv(i,:) = [r', (G*r)']/rho;
end

Jac = inv(Jac_inv);
%% 雅可比矩阵验证
q=@(t)[t;2*t+5;0;0;0;0];
dq = [1;2;0;0;0;0];

Jac = iJacobian(a,b,Len,p,R1);
v1 = Jac*dq;
t = 0:0.1:1;
num = length(t);
p2 = zeros(4,num);
R2 = zeros(3,3,num);
for i =1:num
    [p2(:,i),R2(:,:,i)] = fkine(a,b,Len,q(t(i)),n);
end

%% force
F = [0 0 -1200]';   % 末端载荷
T = [0 100 200]';   % 末端扭矩
L = zeros(3,6);     % 支链向量
Lx = zeros(3,3,6);  % 支链向量叉乘矩阵
bx = zeros(3,3,6);  % 末端铰链叉乘矩阵（用于计算动平台力矩平衡）
Bx = zeros(3,3,6);  % 末端铰链叉乘矩阵（用于计算支链力矩平衡）
m = 0;            % 连杆质量 kg
M = zeros(3,6); M(3,:) = -m*10; % 质量矩阵
N = zeros(3,6);     % 质量引起力矩（参考坐标系下，各支链S1为支点）
for i =1:6
    L(:,i) = p + R1*b(:,i)-a(:,i);
    Lx(:,:,i) = CrossMatrix(L(:,i));
    bx(:,:,i) = CrossMatrix(R1*b(:,i));
    N(:,i) = Lx(:,:,i)/2*M(:,i);
    Bx(:,:,i) = CrossMatrix(p+R1*b(:,i));
end

U1 = eye(18); U2 = eye(18); % 支链力平衡 18*36
U3 = zeros(18);             % 支链力矩平衡 18*18
U4 = blkdiag(Lx(:,:,1),Lx(:,:,2),Lx(:,:,3),Lx(:,:,4),Lx(:,:,5),Lx(:,:,6));  % 支链力矩平衡 18*18
U5 = zeros(6,18);           % 动平台力-力矩平衡（与f1受力无关） 6*18
U6 = [eye(3),    eye(3),    eye(3),    eye(3),    eye(3),    eye(3);        % 动平台力平衡   3*18
    Bx(:,:,1), Bx(:,:,2), Bx(:,:,3), Bx(:,:,4), Bx(:,:,5), Bx(:,:,6) ];   % 动平台力矩平衡 3*18

eq_A = [U1,U2;
    U3,U4;
    U5,U6];

eq_B = [-reshape(M,[18,1]);
    -reshape(N,[18,1]);
    F;
    T+CrossMatrix(p)*F];

x_force = eq_A\eq_B;

drive_force = reshape(x_force(1:18),[3,6]);
plat_force  = reshape(x_force(19:end),[3,6]);

drive_forceref = (Jac'*[F;T])';
check = (drive_force(2,:) - drive_forceref)./drive_force(2,:);

% a_check = zeros(3,6);
for i =1:6
    % a_check = [plat_force(:,i)+drive_force(:,i)+M(:,i)]
    % a_check = N(:,i)+Lx(:,:,i)*plat_force(:,i)
    % a_check = sum(plat_force,2)-F
    % a_check(:,i) = bx(:,:,i)*plat_force(:,i);
end
% aa_check = sum(a_check,2)-T
%% 绘制
vertice = [a b];
vertice(2,7:end)=-Hy*ones(1,6);
face = [1 2 3 4 5 6
    7 8 9 10 11 12];
c = [0;1];
patch("Faces",face(1,:),"Vertices",vertice(:,7:end)','FaceVertexCData',c(2),'FaceColor','flat')
hold on
scatter3([a(1,:),b(1,:)], ...
    [a(2,:),-Hy*ones(1,6)], ...
    [a(3,:),b(3,:)], ...
    70,'red','filled')
for i = 1:6
    plot3([a(1,i),b(1,i)],[a(2,i),-Hy],[a(3,i),b(3,i)],'k','LineWidth',1)
    plot3([a(1,i),a(1,i)],[-Hy,Hy],[a(3,i),a(3,i)],':b','LineWidth',2)
end
xlabel('X轴')
ylabel('Y轴')
zlabel('Z轴')

hold off
view(3)
%% function addition
function [X]=initialization(pop,ub,lb,dim)
    for i=1:pop
        for j=1:dim
            X(i,j)=(ub(j)-lb(j))*rand()+lb(j);%在限定的  
        end 
    end
end
 
function fitness=fun(X)
    % Hy = 150;              % Initial height
    n = 1;
    A_L = 400;
    A_d = X(1);
    B_L = X(2);
    B_d = X(3);
    Hy  = X(4);
    p0 = [0,-Hy,0]'/n;
    R0 = eye(3);
    % fitness = 0;
    [a,b]=Initial_Config(A_d, A_L, B_d, B_L,0,n);

    Len = vecnorm(repmat(p0,1,6) + R0*b - a);
    
    % Len = zeros(1,6);
    % for i = 1:6
    %     Len(i) = norm(R0*b(:,i)+p0-a(:,i)); % 支链长度，单位与n保持一致
    % end

    % 定义参数范围
    x_values = -25:5:25;
    y_values = -25:5:25;
    z_values = -25:5:25;  
    ry_values = -3:2:3;
    rx_values = -3:2:3;
    rz_values = -30:10:30;

    % 生成所有参数组合（使用ndgrid展开）
    [X, Y, Z, RY, RZ, RX] = ndgrid(x_values, y_values, z_values, ry_values, rz_values, rx_values);

    % 转换为列向量
    params = [X(:), Y(:), Z(:), RY(:), RZ(:), RX(:)];

    % 初始化存储结果的数组
    num_params = numel(params)/6;  % 总参数组合数
    fitness_values = zeros(num_params, 1);

    % 启动并行池（如果尚未启动）
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    sdot = @(s1,s2)sum(s1(1:3,:).*s2(4:6,:)+s1(4:6,:).*s2(1:3,:),1);
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
        sk = ikine(a,b,Len,pk,Rk)
        ak = a;
        ak(2,:) = sk;
        Lk = repmat(pk,1,6) + Rk*b - ak;
        
        C_phik = ([0,-1,0]*Lk)./vecnorm(Lk,2,1);
        lambda_k = abs(C_phik);

        % 计算雅可比矩阵及其条件数
        [Jac,~] = iJacobian(a, b, Len, pk, Rk); % 计算雅可比矩阵，输入单位均一致。
        SO_k = Jac + [skew(pk)*Jac(4:6,:);zeros(3,6)];
        SO_k = SO_k./vecnorm(SO_k(4:6,:),2,1);
        ST_k = [zeros(3,6);Lk];

        ita_k = abs(sdot(SO_k,ST_k))./ (vecnorm(SO_k(1:3,:),2,1).*vecnorm(ST_k(4:6,:),2,1));

        
        fitness_values(i) = min([ita_k,lambda_k]);
    end

    % 获取全局最优解
    [fitness, ~] = max(1./fitness_values);
    if isempty(fitness)
        disp('something wrong')
        fitness = 0;
    end
end

 
function [X]=BoundaryCheck(X,ub,lb,dim)
    for i=1:dim
        if X(i)>ub(i)
            X(i)=ub(i);
        end
        if X(i)<lb(i)
            X(i)=lb(i);
        end
    end
end
 
function [Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,fobj,vmax,vmin,maxIter)
    c1=2.0;
    c2=2.0;
    w_max = 2;
    w_min = 0.4;
    V=initialization(pop,vmax,vmin,dim);
    X=initialization(pop,ub,lb,dim);
    fitness=zeros(1,pop);

    for i=1:pop
        fitness(i)=fobj(X(i,:));
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
            V(i,:)=BoundaryCheck(V(i,:),vmax,vmin,dim);
            Xnew(i,:)=X(i,:)+V(i,:);
            Xnew(i,:)=BoundaryCheck(Xnew(i,:),ub,lb,dim); 

            fitnessNew(i)=fobj(Xnew(i,:));
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