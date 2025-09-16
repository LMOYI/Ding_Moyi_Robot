classdef Stewart
    %Stewart This file is a module for describing the kinematics of parallel robot Stewart
    %   Author: Jinke Li.(Based on Weijia Zhang's version)
    %   version 1.0
    %   using SE3 1.0
    
    properties (Constant)
        myse3=mySE3();      % SE3 module
        myudq = myUDQ();
    end
    properties
        n;              % int: n=1 in mm, n=1e-3, in m
        axes;           % cell<1,3> -> cell<1,*> -> matrix<6,1> : axes represented in the base frame
        ra;             % double : the radius of base platform
        rb              % double : ** of moving platform
        L;              % double
        H;              % double
        g0;             % matrix<4,4> : initialPose of the moving platform
        q0;             % matrix<3,1> : initial variables of active joints
        jointType;      % matrix<1,36> -> char : axes types
        joint2axes;     % cell(1,36) -> cell(1,*) -> matrix<1,2> : the mappings from joints to their corresponding axes
        axes2joint;     % cell(1,num_of_chains) -> matrix<1,num_of_axes> : the sequence for active and passive joints seperately
        jointPoint;     % cell<1,6> -> matrix<6,1> : specified points on the axes
        jointBound;     % cell(1,6) -> matrix<1,2> : bounds for joint variables (only bounds for active joints matter in kinematic calibration)
        jointActive;    % matrix<1,6> : the sequence of active joints
        jointPassive;   % matrix<1,30> : the sequence of passive joints
        axesActive;     % cell(1,num_of_chains) -> matrix<1,num_of_active_axes>
        axesPassive;    % cell(1,num_of_chains) -> matrix<1,num_of_passive_axes>
        A;              % cell<1,6> -> matrix<3,1> : intersection points of the base platform and the active joint axes
        B;              % cell<1,6> -> matrix<3,1> : the joint centers in the base frame at the initial pose
        C;              % cell<1,6> -> matrix<3,1>
    end
    
    
    
    methods
        function [obj] = Stewart(ra,rb,L,n)
            %Stewart 构造此类的实例
            %   input:
            %   angle_bound::double
            %   position_bound::double
            %   n::int, n==1 unit mm, n==1e-3, unit m
            %   output:
            %   obj::class(PSR2PRU)
            %   sa::cell<1,3> -> cell<1,*> -> matrix<6,1>
            %   ga0::matrix<4,4>
            mySE3=obj.myse3;
            % r= 190*n;
            % L= 474*n;

            theta=asin(55/190);%基座平台上连接点相对于中心的角度，这里取用左边第一个连接点
            % H=sqrt(L^2-4*r^2*sin(pi/6-theta)^2);%基座平台与移动平台之间的初始高度
            
            A1=ra*[cos(pi/3-theta), sin(pi/3-theta), 0]';
            A2=ra*[cos(pi/3+theta), sin(pi/3+theta), 0]';
            A3=ra*[cos(pi-theta), sin(pi-theta), 0]';
            A4=ra*[cos(pi+theta), sin(pi+theta), 0]';
            A5=ra*[cos(-pi/3-theta), sin(-pi/3-theta), 0]';
            A6=ra*[cos(-pi/3+theta), sin(-pi/3+theta), 0]';
            A={A1,A2,A3,A4,A5,A6};%基座六个支点 x轴水平向右，y轴水平向上，四个区间逆时针定位
            
            C1=rb*[cos(theta), sin(theta), 0]';
            C2=rb*[cos(2*pi/3-theta), sin(2*pi/3-theta), 0]';
            C3=rb*[cos(2*pi/3+theta), sin(2*pi/3+theta), 0]';
            C4=rb*[cos(-2*pi/3-theta), sin(-2*pi/3-theta), 0]';
            C5=rb*[cos(-2*pi/3+theta), sin(-2*pi/3+theta), 0]';
            C6=rb*[cos(-theta), sin(-theta), 0]';
            C={C1,C2,C3,C4,C5,C6};%动平台坐标系在基座坐标系的映射位置
            
            Projected_length = norm(C1-A1);
            H = sqrt(L^2-Projected_length^2);
            
            B1=rb*[cos(theta), sin(theta), 0]'+[0, 0, H]';
            B2=rb*[cos(2*pi/3-theta), sin(2*pi/3-theta), 0]'+[0, 0, H]';
            B3=rb*[cos(2*pi/3+theta), sin(2*pi/3+theta), 0]'+[0, 0, H]';
            B4=rb*[cos(-2*pi/3-theta), sin(-2*pi/3-theta), 0]'+[0, 0, H]';
            B5=rb*[cos(-2*pi/3+theta), sin(-2*pi/3+theta), 0]'+[0, 0, H]';
            B6=rb*[cos(-theta), sin(-theta), 0]'+[0, 0, H]';
            B={B1,B2,B3,B4,B5,B6};%动平台高度，x轴水平向右，y轴水平向上，四个区间逆时针定位
            

            
            %  initial pose末端动平台
            Ro=eye(3);
            po=[0,0,H]';
            g0=[Ro po; 0 0 0 1];
            

            % chains
            sn=cell(1,6);%sn 是一个 1x6 的单元数组，表示 Stewart 平台的 6 个支腿。
            for c=1:6
                sn{c}=cell(1,6);%对于每个支腿 c，sn{c} 是一个 1x6 的单元数组，表示该支腿的 6 个关节轴。
                so1=[1,0,0]';
                so2=[0,1,0]';
                so3=(B{c}-A{c})/norm(B{c}-A{c});
                so4=[1,0,0]';
                so5=[0,1,0]';
                so6=[0,0,1]';
                
                sn{c}{1}=[so1;cross(A{c},so1)];
                sn{c}{2}=[so2;cross(A{c},so2)];
                sn{c}{3}=[zeros(3,1);so3];
                sn{c}{4}=[so4;cross(B{c},so4)];
                sn{c}{5}=[so5;cross(B{c},so5)];
                sn{c}{6}=[so6;cross(B{c},so6)];%旋量形式
            end
            % arbitrarily set the points on the joints at the initial pose
            %第一个支腿的前三个关节点位于平台A，后三个关节点位于平台B   不考虑中间平台
            pn={A1,A1,A1,B1,B1,B1,A2,A2,A2,B2,B2,B2,A3,A3,A3,B3,B3,B3,A4,A4,A4,B4,B4,B4,A5,A5,A5,B5,B5,B5,A6,A6,A6,B6,B6,B6};
            % set the motion ranges of active joints%移动关节运动范围
            d1_bound=[0, 150]*n;
            d2_bound=[0, 150]*n;
            d3_bound=[0, 150]*n;
            d4_bound=[0, 150]*n;
            d5_bound=[0, 150]*n;
            d6_bound=[0, 150]*n;
            jBound={d1_bound,d2_bound,d3_bound,d4_bound,d5_bound,d6_bound};
            % joint types
            jType=zeros(1,36);
            % joint control modes
            jActive=zeros(1,6);%6个主动关节，30个被动关节
            jPassive=zeros(1,30);
            % mapping from joints to axes
            joint2axes=cell(1,36);%一种索引，对应第几个支腿的第几个轴
            for c=1:6
                jType(6*(c-1)+1:6*c)=['r','r','p','r','r','r'];%每个支腿的类型
                jActive(c)=6*(c-1)+3;%主动腿的索引
                jPassive(5*(c-1)+1:5*c)=[6*(c-1)+1,6*(c-1)+2,6*(c-1)+4,6*(c-1)+5,6*(c-1)+6];%被动腿的索引
                for a=1:6
                    joint2axes{6*(c-1)+a}={[c,a]};
                end
            end
            axes2joint=cell(1,6);
            axes2joint{1}=[1,2,1,3,4,5];
            axes2joint{2}=[6,7,2,8,9,10];
            axes2joint{3}=[11,12,3,13,14,15];
            axes2joint{4}=[16,17,4,18,19,20];
            axes2joint{5}=[21,22,5,23,24,25];
            axes2joint{6}=[26,27,6,28,29,30];
            axesActive=cell(1,6);
            axesActive{1}=[3];
            axesActive{2}=[3];
            axesActive{3}=[3];
            axesActive{4}=[3];
            axesActive{5}=[3];
            axesActive{6}=[3];
            axesPassive=cell(1,6);
            axesPassive{1}=[1,2,4,5,6];
            axesPassive{2}=[1,2,4,5,6];
            axesPassive{3}=[1,2,4,5,6];
            axesPassive{4}=[1,2,4,5,6];
            axesPassive{5}=[1,2,4,5,6];
            axesPassive{6}=[1,2,4,5,6];
            % model parameters
            obj.axes=sn;
            obj.ra=ra;
            obj.rb=rb;
            obj.L=L;
            obj.H=H;
            obj.g0=g0;
            obj.q0=L*ones(6,1);
            obj.jointType=jType;
            obj.joint2axes=joint2axes;
            obj.axes2joint=axes2joint;
            obj.jointPoint=pn;
            obj.jointBound=jBound;
            obj.jointActive=jActive;
            obj.jointPassive=jPassive;
            obj.axesActive=axesActive;
            obj.axesPassive=axesPassive;
            obj.A=A;
            obj.B=B;
            obj.C=C;
            obj.n=n;
        end
        
        function [q] = ikine(obj,g)% 简化逆运动学
            %ikines 简化的反向运动学
            %   input:
            %   g::matrix<6,6>
            %   output:
            %   q::matrix<6,1>
            q=zeros(6,1);
            for i=1:6
                P=g(1:3,1:3)*obj.C{i}+g(1:3,4);
                q(i)=norm(P-obj.A{i});
            end
            q=q-obj.q0;
        end
        
        function [Q,residuals] = ikines(obj,ga)
            %ikine 反向运动学
            %   input:
            %   g::matrix<6,6> 目标输入角度
            %   output:
            %   Q::cell(1,6)
            mySE3=obj.myse3();
            Q=cell(1,6);%Q为输出的关节角度
            residuals=zeros(6,1);%residual为残差，根据这个残差进行迭代
            axes=obj.axes;%表示每条支链每个关节的旋量
            for c=1:6%c代表每条腿
                Q{c}=[0,0,0,0,0,0]';    %每条支链上各关节移动距离
            end
            Lmax=0.1;%把每次迭代求出的dq模长控制在一定范围，即每次关节角度改变控制在一定范围内
            semax=0.02;%用于求解步长的一个参数
            Te=mySE3.Tinv(obj.g0)*ga;%inv(g0)*ga 表示一个误差变换SE(3)，其中误差源于初始位姿和目标位姿的差距
            se=mySE3.log2s(Te);%转换为误差旋量
            steps=ceil(max(abs(se)/semax));%对误差进行分步
            seL=se/steps;
            for t=1:steps
                gtmp=obj.g0*mySE3.exp2T(seL*t);%每一步初始位姿和增量组合成新目标位姿
                for c=1:6 %此下都是对每个腿分析
                    g=cell(1,6);%单关节的变换矩阵
                    G=cell(1,7);%每个关节关于首个关节的变换关系
                    for k=1:10
                        for a=1:6 %a代表每个等效关节
                            g{a}=mySE3.exp2T(axes{c}{a}*Q{c}(a));%axes是每个关节的旋量 Q是旋量旋转的角度
                        end
                        G{1}=eye(4);
                        for a=2:7
                            G{a}=G{a-1}*g{a-1};%表示支链的变换关系
                        end
                        Gst=G{end}*obj.g0;%带初始位置的本支链腿的当前位姿
                        y=[eye(3),zeros(3);-mySE3.Sw(Gst(1:3,4)),eye(3)]*mySE3.log2s(gtmp*mySE3.Tinv(Gst));
                        %y是末端误差向量，不直接调用正向运动学避免迭代，把gtmp*SE3.Tinv(Gst)（误差SE3）转为旋量，然后通过伴随变换映射到末端
                        %旋转的部分没有变化，平移的部分要减去旋转对平移的影响，v-pxv
                        residuals(c)=norm(y);%residual是残差，表示当前位姿与目标位姿的差异
                        Jq=zeros(6,1);%雅可比矩阵，对每个关节旋量运用伴随变换映射到基座标系，然后伴随变换映射到末端
                        for a=1:6
                            Jq(:,a)=[eye(3),zeros(3);-mySE3.Sw(Gst(1:3,4)),eye(3)]*mySE3.Ad(G{a})*axes{c}{a};
                        end
                        dq=Jq\y;%x_dot=Jq_dot
                        if max(abs(dq))>Lmax
                            dq=dq/Lmax;
                        end
                        Q{c}=Q{c}+dq;
                        if norm(y)<1e-10
                            break;
                        end
                    end
                end
            end
        end
        
        function [g] = fkine(obj,q,x0)
            %fkine 正向运动学
            %   input:
            %   q::matrix<1,6>
            %   output:
            %   g::matrix<6,6>
            udqtool=obj.myudq;
            a = [obj.A{1},obj.A{2},obj.A{3},obj.A{4},obj.A{5},obj.A{6}];
            b = [obj.C{1},obj.C{2},obj.C{3},obj.C{4},obj.C{5},obj.C{6}];
            
            % 对偶四元数方法
            % 此时 x = [xi;lambda]
            % xi = [x1,x2,x3,x0];
            % [w,theta] = Quaternion(R);
            % Representation of pose of mobile platform

            conj = @(x)[-x(1:3);x(4)];
            skew = @(x)[0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
            
            % 原理公式
            % E_1 = p + Rb =  d + l_i;
            % E_2(theta)*xi = lambda + xi*ai + ci*xi
            % norm(E_2)^2 = **** = norm(l_i)^2 = x'Qx/2 - C_i
            q = reshape(q,1,6);
            Len = obj.L + q;
            Q = zeros(8,8,8);
            Q(1:4, 1:4, 7)=2*eye(4); % xi 模长为1
            Q(1:4, 5:8, 8)=eye(4);
            Q(5:8, 1:4, 8)=eye(4);   % xi 和 lambda 正交

            q_C = [Len.^2,1,0]';
            %　a-base, b-mobile platform
            for i = 1:6

                % Qi11-- 4*4-- 0.5*xi'Qxi
                b1 = [-skew(b(:,i)), b(:,i); -b(:,i)', 0];
                c1 = [skew(-a(:,i)),-a(:,i);  a(:,i)', 0];
                Qi11 = 2*(b1'*b1) + 2*(c1'*c1) + 4*(c1'*b1);
                Qi21 = 2*b1 + 2*c1;
                Qi12 = Qi21';
                Qi22 = 2*eye(4);
                Q(:,:,i) = [Qi11,Qi12;Qi21,Qi22];
            end

            flag = 1;
            if nargin == 2 || isempty(x0)
                x0 = [0,0,0,1,0,0,1,0]';
            end
            xi_0 = x0(1:4);           % 四元数-旋转相关
            lambda_0 = x0(5:8);      % 四元数-平移相关
            X = [xi_0;lambda_0];

            while norm(flag)>1e-10
                J = [Q(:,:,1)*X, Q(:,:,2)*X, Q(:,:,3)*X, Q(:,:,4)*X, Q(:,:,5)*X, Q(:,:,6)*X, Q(:,:,7)*X, Q(:,:,8)*X]';
                % % delta_X = (J'*J)\(J'*C);
                delta_X = pinv(J)*q_C;
                % delta_X = J\C;
                X1 = 0.5*X + delta_X;
                flag =  X1-X;
                X=X1;
                % norm(flag)
            end

            % X0 = reshape(X,8,1);
            xi  = X(1:4);
            lambda = X(5:8);

            % xi = [w.*sind(theta/2);cosd(theta/2)];
            % lambda := p*xi
            theta = 2*acos(xi(4));
            if theta~=0
                w = xi(1:3)/sin(theta/2);
            else
                w = [0;1;0];
            end
            
            % 由于该方法中四元数定义和udqtool中定义不一致，所以需要做如下转化
            lambda = lambda([4,1,2,3]);
            xi_conj = conj(xi);
            xi_conj = xi_conj([4,1,2,3]);
            % % 保持输出一致[整个过程按m计算]
            % if obj.n == 1e-3   % 输出 m，
            %     p = udqtool.Q_Times(lambda,xi_conj)*1;
            % else        % 输出 mm
            %     p = udqtool.Q_Times(lambda,xi_conj);
            % end
            p = udqtool.Q_Times(lambda,xi_conj);
            % R = rotation_matrix(w,theta);
            R = skew(w)*sin(theta) + skew(w)*skew(w)*(1-cos(theta)) + eye(3);
            g = [R,p(2:4);0,0,0,1];

        end
        
        function J_inv = Jacobian(obj,g)
            % J dx = dq, dx = [V;W]
            % J_inv 是几何雅可比 J_inv * dq = dx = [v;w]
            % qF' * dq = T' * dx = T' * J_inv*dq
            % qF' = T'*J_inv
            R = g(1:3,1:3);
            unit = obj.n*1e3; % 统一单位为m
            
            p = g(1:3,4)/unit;
            b = [obj.C{1},obj.C{2},obj.C{3},obj.C{4},obj.C{5},obj.C{6}]/unit;
            a = [obj.A{1},obj.A{2},obj.A{3},obj.A{4},obj.A{5},obj.A{6}]/unit;
            Link_Vec = R*b+repmat(p,1,6) - a;
            L_unit = Link_Vec*diag(1./vecnorm(Link_Vec));
            J = [L_unit', cross(b'*R',L_unit',2)];
            J_inv = inv(J);
        end

        function index = Index(q,type)
            % type = 1, index = LTI
            % type = 0, index = 
            
        end

        function qF = force(obj,g,F,T)
            unit = obj.n*1e3; % 统一单位为m
            R = g(1:3,1:3);
            p = g(1:3,4)/unit;
            
            a = [obj.A{1},obj.A{2},obj.A{3},obj.A{4},obj.A{5},obj.A{6}]/unit;
            b = [obj.C{1},obj.C{2},obj.C{3},obj.C{4},obj.C{5},obj.C{6}]/unit;
            Link_Vec = (R*b+p-a);
            L_unit = Link_Vec*diag(1./vecnorm(Link_Vec));
            % Sum qFi * L_unit_i = F           
                % L_unit*qF = F
            % Sum cross(ai, qFi*Link_Veci) = T + cross(p,F)   
                % a_hat * Link_Vec * qF = T + p_hat*F
            Hq = [L_unit;cross(a,L_unit)];
            qF = Hq\[F;T+cross(p,F)];
        
        end

        function Plot_Robot(obj,g)
            R = g(1:3,1:3);
            p = g(1:3,4);
            a = [obj.A{1},obj.A{2},obj.A{3},obj.A{4},obj.A{5},obj.A{6}];
            b = [obj.C{1},obj.C{2},obj.C{3},obj.C{4},obj.C{5},obj.C{6}];

            % 计算连杆向量 
            Link_Vec = (R*b + p - a);
            
            % 这是绘制动平台本身所必需的
            b_world = R * b + p;
            
            h_out = []; % 初始化输出句柄数组
            
            % a. 绘制基座平台 (半透明灰色)
            h = patch('Vertices', a', 'Faces', 1:6, ...
                      'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
            hold on
            
            h_out = [h_out; h]; % 将句柄存入数组
            
            % b. 绘制动平台 (半透明蓝色)
            h = patch('Vertices', b_world', 'Faces', 1:6, ...
                      'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.6);
            h_out = [h_out; h]; % 将句柄存入数组
            
            % c. 绘制6个连杆/支腿 (深蓝色实线)
            for i = 1:6
                % 每条连杆的坐标是 [基座连接点, 动平台世界坐标连接点]
                leg_coords = [a(:, i), b_world(:, i)];
                
                h = plot3(leg_coords(1,:), leg_coords(2,:), leg_coords(3,:), ...
                          'Color', [0.07, 0.16, 0.44], 'LineWidth', 4);
                h_out = [h_out; h]; % 将句柄存入数组
            end

        end
        
        function Robot_Check(obj,Lim_Kine,type,Force)

            if isempty(Force) || nargin == 3
                F = [0;0;1200];
                T = [0;0;0];
            else
                F = Force(1:3);
                T = Force(4:6);
            end
            sdot = @(s1,s2)sum(s1(1:3,:).*s2(4:6,:)+s1(4:6,:).*s2(1:3,:),1); % 定义互易积
            xlim = Lim_Kine(:,1);
            ylim = Lim_Kine(:,2);
            zlim = Lim_Kine(:,3);
            rxlim = Lim_Kine(:,4);
            rylim = Lim_Kine(:,5);
            rzlim = Lim_Kine(:,6);
            %定义参数范围
            x_values = xlim(1):5:xlim(2);
            y_values = ylim(1):5:ylim(2);
            z_values = zlim(1):5:zlim(2);  % 注意：这个范围实际只会生成25
            ry_values = rylim(1):5:rylim(2);
            rx_values = rxlim(1):5:rxlim(2);
            rz_values = rzlim(1):5:rzlim(2);
            nx= length(x_values);
            ny= length(y_values);
            nz= length(z_values);
            nrx= length(rx_values);
            nry= length(ry_values);
            nrz= length(rz_values);
            
            % 生成所有参数组合（使用ndgrid展开）
            [X, Y, Z, RY, RZ, RX] = ndgrid(x_values, y_values, z_values, ry_values, rz_values, rx_values);
            
            % 转换为列向量
            params = [X(:), Y(:), Z(:), RY(:), RZ(:), RX(:)];
            
            %初始化存储结果的数组
            num_params = numel(params)/6;  % 总参数组合数
            fitness_values1 = zeros(num_params, 1);% 条件数
            fitness_values2 = zeros(num_params, 1);% LTI
            F_joint = zeros(6,num_params);
            F_joint_ref = zeros(6,num_params);
            S_length = zeros(num_params, 6);
            %启动并行池（如果尚未启动）
            if isempty(gcp('nocreate'))
                parpool('local');
            end
            
            a = [obj.A{1},obj.A{2},obj.A{3},obj.A{4},obj.A{5},obj.A{6}];
            b = [obj.C{1},obj.C{2},obj.C{3},obj.C{4},obj.C{5},obj.C{6}];
            %并行计算
            parfor i = 1:num_params
                % 解包参数
                x = params(i,1);
                y = params(i,2);
                z = params(i,3);
                ry = params(i,4);
                rz = params(i,5);
                rx = params(i,6);
            
                % 计算变换矩阵
                pk = [x; y; z];
                Rk = rotz(rz, "deg") * roty(ry, "deg") * rotx(rx, "deg");
                gk = [Rk,pk;0,0,0,1];
                s = obj.ikine(gk);
                S_length(i,:) = s + obj.L;
                % 计算雅可比矩阵及其条件数
                
                % force_x(i,:) = f(1,:);
                % force_y(i,:) = f(2,:);
                % force_z(i,:) = f(3,:);
                Jac = obj.Jacobian(gk);
                fitness_values1(i) = cond(Jac);

                Lk = repmat(pk,1,6) + Rk*b - a;

                % 输入传递指标 lambda
                % TWS（传递力旋量） 与 输入旋量 互易积
                % TWS 传递到末端平台的力
                % C_phik = (Lk*Lk)./vecnorm(Lk,2,1);
                % lambda_k = abs(C_phik);
                lambda_k = 1;
                % 计算几何雅可比矩阵及其条件数
                Jac = obj.Jacobian(gk);
                % fitness_values(i) = cond(Jac);

                % 输出传递指标
                % TWS 和 输出旋量 互易积
                % Jac*q=v = [v;w]
                SO_k = zeros(6,6);
                for j = 1:6
                    Tk = obj.e2g(Jac(:,j),1);
                    dTk_hat = logm(Tk);
                    vk = dTk_hat(1:3,4);
                    wk = [dTk_hat(3,2);dTk_hat(1,3);dTk_hat(2,1)];
                    SO_k(:,j) = [vk;wk];
                end
                % SO_k = Jac + [skew(pk)*Jac(4:6,:);zeros(3,6)]; 
                SO_k = SO_k./vecnorm(SO_k(4:6,:),2,1);

                ST_k = [zeros(3,6);Lk]./vecnorm(Lk,2,1);
                ita_k = abs(sdot(SO_k,ST_k))./ (vecnorm(SO_k(1:3,:),2,1).*vecnorm(ST_k(4:6,:),2,1));

                fitness_values2(i) = min([ita_k,lambda_k]);
                % 计算力
                qf = obj.force(gk,F,T);
                F_joint(:,i) = qf;
                F_joint_ref(:,i) = Jac'*[F;T];
            end

            %获取全局最优解
            [CondJ, ~] = max(fitness_values1);
            N_all = reshape(fitness_values1,nx*ny*nz,nrx*nry*nrz);
            N2_all = reshape(fitness_values2,nx*ny*nz,nrx*nry*nrz);
            N = max(N_all,[],2);
            N2 = max(N2_all,[],2);

            [X,Y,Z] = ndgrid(x_values, y_values, z_values);
            N_matrix= reshape(N,nx,ny,nz);
            N2_matrix= reshape(N2,nx,ny,nz);

            % --- 可视化
            if type ==1
                
                figure('Name', '方法A: 三维散点图', 'Position', [100, 100, 800, 600]);
                scatter3(X(:), Y(:), Z(:), 30, N_matrix(:), 'filled');
                
                title('工作空间条件数分布 - 散点图');
                xlabel('X (mm)');
                ylabel('Y (mm)');
                zlabel('Z (mm)');
                h = colorbar; % 添加颜色条
                ylabel(h, 'Condition Number'); % 给颜色条添加标签
                clim([min(N) max(N)]); % 统一颜色范围
                view(3); % 设置为三维视角
                axis equal;
                grid on;

                figure('Name', '方法A: 三维散点图', 'Position', [100, 100, 800, 600]);
                scatter3(X(:), Y(:), Z(:), 30, N2_matrix(:), 'filled');
                
                title('工作空间LTI分布 - 散点图');
                xlabel('X (mm)');
                ylabel('Y (mm)');
                zlabel('Z (mm)');
                h = colorbar; % 添加颜色条
                ylabel(h, 'Condition Number'); % 给颜色条添加标签
                clim([min(N2) max(N2)]); % 统一颜色范围
                view(3); % 设置为三维视角
                axis equal;
                grid on;
            elseif type == 2
                
                figure('Name', '方法B: 切片图', 'Position', [150, 150, 800, 600]);
                % 定义要显示切片的位置
                slice_x = [xlim(1), (xlim(1)+xlim(2))/2, xlim(2)]; % X方向切三刀
                slice_y = [ylim(1), ylim(2)];                % Y方向切两刀
                slice_z = [zlim(1), (zlim(1)+zlim(2))/2];       % Z方向切两刀

                % 注意：slice函数需要meshgrid格式的网格，我们这里转换一下
                [X_mesh, Y_mesh, Z_mesh] = meshgrid(x_values, y_values, z_values);
                % 相应地，数据矩阵也需要调整维度顺序
                N_matrix_mesh = permute(N_matrix, [2 1 3]);

                slice(X_mesh, Y_mesh, Z_mesh, N_matrix_mesh, slice_x, slice_y, slice_z);

                title('工作空间条件数分布 - 切片图');
                xlabel('X (mm)');
                ylabel('Y (mm)');
                zlabel('Z (mm)');
                h = colorbar;
                ylabel(h, 'Condition Number');
                clim([min(N) max(N)]);
                view(-35, 25); % 调整视角
                axis equal;
                grid on;

            % elseif type == 3
            %     figure('Name', '方法C: 等值面图', 'Position', [200, 200, 800, 600]);
            % 
            %     % 定义您关心的条件数值（等值面阈值）
            %     cond_threshold_low = 14.4;
            %     cond_threshold_high = 15.1;
            % 
            %     % 绘制条件数等于 low 值的等值面
            %     p1 = patch(isosurface(X, Y, Z, N_matrix, cond_threshold_low));
            %     isonormals(X, Y, Z, N_matrix, p1); % 计算法线，让光照更真实
            %     p1.FaceColor = 'green';
            %     p1.EdgeColor = 'none';
            %     p1.FaceAlpha = 0.5; % 设置透明度
            % 
            %     hold on; % 在同一张图上继续画
            % 
            %     % 绘制条件数等于 high 值的等值面
            %     p2 = patch(isosurface(X, Y, Z, N_matrix, cond_threshold_high));
            %     isonormals(X, Y, Z, N_matrix, p2);
            %     p2.FaceColor = 'red';
            %     p2.EdgeColor = 'none';
            %     p2.FaceAlpha = 0.5;
            % 
            %     % --- 图形美化 ---
            %     title('工作空间条件数分布 - 等值面');
            %     xlabel('X (mm)');
            %     ylabel('Y (mm)');
            %     zlabel('Z (mm)');
            %     legend([p1, p2], {['Cond = ', num2str(cond_threshold_low)], ['Cond = ', num2str(cond_threshold_high)]});
            %     grid on;
            %     axis equal;
            %     view(3);
            %     camlight; % 添加光源
            %     lighting gouraud; % 设置光照模式
            end


            % contour(xlim(1):2:xlim(2),ylim(1):2:ylim(2),N_matrix(:,:,5)');
            % xlabel("x/mm");
            % ylabel("y/mm");
            % title("最大雅可比条件数图");

            F1 = reshape(F_joint(1,:),1,num_params);
            F2 = reshape(F_joint(2,:),1,num_params);
            F3 = reshape(F_joint(3,:),1,num_params);
            F4 = reshape(F_joint(4,:),1,num_params);
            F5 = reshape(F_joint(5,:),1,num_params);
            F6 = reshape(F_joint(6,:),1,num_params);
            
            disp('关节1最大受力')
            disp(max((F1)))
            disp('关节2最大受力')
            disp(max((F2)))
            disp('关节3最大受力')
            disp(max((F3)))
            disp('关节4最大受力')
            disp(max((F4)))
            disp('关节5最大受力')
            disp(max((F5)))
            disp('关节6最大受力')
            disp(max((F6)))
            
            disp('关节1最大驱动力')
            disp(max(abs(F_joint_ref(1,:))))
            disp('关节2最大驱动力')
            disp(max(abs(F_joint_ref(2,:))))
            disp('关节3最大驱动力')
            disp(max(abs(F_joint_ref(3,:))))
            disp('关节4最大驱动力')
            disp(max(abs(F_joint_ref(4,:))))
            disp('关节5最大驱动力')
            disp(max(abs(F_joint_ref(5,:))))
            disp('关节6最大驱动力')
            disp(max(abs(F_joint_ref(6,:))))
            
            S = (max(S_length,[],1) - min(S_length,[],1))*obj.n
            S_1 = max(S_length,[],1)*obj.n
            S_2 = min(S_length,[],1)*obj.n
        end

        function WorkSpace_Plot(type)

            
        end

        function Robot_JointsRange(obj,Lim_Kine)
            % 用于计算输入范围内各关节角度范围
            % Lim_Kine -> matrix<2,6>  -> [xlim,ylim,zlim,rxlim,rylim,rzlim]

            xlim = Lim_Kine(:,1);
            ylim = Lim_Kine(:,2);
            zlim = Lim_Kine(:,3);
            rxlim = Lim_Kine(:,4);
            rylim = Lim_Kine(:,5);
            rzlim = Lim_Kine(:,6);
            %定义参数范围
            x_values = xlim(1):2:xlim(2);
            y_values = ylim(1):2:ylim(2);
            z_values = zlim(1):2:zlim(2); 
            ry_values = rylim(1):2:rylim(2);
            rx_values = rxlim(1):2:rxlim(2);
            rz_values = rzlim(1):2:rzlim(2);
            nx= length(x_values);
            ny= length(y_values);
            nz= length(z_values);
            nrx= length(rx_values);
            nry= length(ry_values);
            nrz= length(rz_values);
            
            % 生成所有参数组合（使用ndgrid展开）
            [X, Y, Z, RY, RZ, RX] = ndgrid(x_values, y_values, z_values, ry_values, rz_values, rx_values);
            
            % 转换为列向量
            params = [X(:), Y(:), Z(:), RY(:), RZ(:), RX(:)];
            
            %初始化存储结果的数组
            num_params = numel(params)/6;  % 总参数组合数
            
            %启动并行池（如果尚未启动）
            if isempty(gcp('nocreate'))
                parpool('local');
            end
            
            a = [obj.A{1},obj.A{2},obj.A{3},obj.A{4},obj.A{5},obj.A{6}];
            b = [obj.C{1},obj.C{2},obj.C{3},obj.C{4},obj.C{5},obj.C{6}];
            %并行计算
            parfor i = 1:num_params
                % 解包参数
                x = params(i,1);
                y = params(i,2);
                z = params(i,3);
                ry = params(i,4);
                rz = params(i,5);
                rx = params(i,6);
            
                % 计算变换矩阵
                pk = [x; y; z];
                Rk = rotz(rz, "deg") * roty(ry, "deg") * rotx(rx, "deg");
                gk = [Rk,pk;0,0,0,1];
                s = obj.ikines(gk);
                
                
            end
        end

            function g = e2g(obj,e,type)
                if type == 1 % 广义旋转
                    p = zeros(3,1);
                    R = eye(3);
                    p(:) = e(1:3);
                    w_hat = [0,-e(6),e(5); e(6),0,-e(4);-e(5),e(4),0];
                    R(:,:) = expm(w_hat);
                    g = [R,p;0,0,0,1];
                end
            end
    end


end