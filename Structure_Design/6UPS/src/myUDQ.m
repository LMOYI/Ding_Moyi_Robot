classdef myUDQ

    properties
        % p = [p0,p1,p2,p3]' p0是实部
        % P = [p1;p2];
    end
    methods
        function obj = myUDQ()
            % obj.P = zeros(4,n);
            % obj.D = zeros(4,n);
            % obj.H = [obj.P;obj.D];
            %myUDQ 构造此类的实例
        end

        function [b,c] = PD(obj,a)
            %PD 计算a的主部（Primal）和对偶部（Dual）
            %Input:     a -> matrix<8,Num> -> UDQ
            %Output:    b -> matric<4,Num> -> Primal part of a
            %           c -> matrix<4,Num> -> Dual part of a
            b = a(1:4,:);
            c = a(5:8,:);
        end

        function h = UDQ_Exp(obj,s)
            %UDQ_Exp 计算s的指数映射，切平面映射到流形上？得到UDQ形式
            %%%Input    s -> matrix<6,n> -> s = phi*n/2 + epslon*p/2 -> DQ纯虚部
            %%%Output:  h -> matrix<8,n> -> h = r + epslon*p*r/2 -> UDQ
            n = size(s,2);
            theta = zeros(1,n);
            theta_s = zeros(1,n);
            for i = 1:n
                if norm(s(1:3,i))>1e-10
                    theta(i) = norm(s(1:3,i));
                    theta_s(i) = sin(theta(i))./theta(i);
                else
                    theta(i) = 0;
                    theta_s(i) = 1;
                end
            end
            % h = zeros(8,n);
            h = [cos(theta);
                s(1,:).*theta_s;
                s(2,:).*theta_s;
                s(3,:).*theta_s;
               -s(1,:).*s(4,:).*theta_s - s(2,:).*s(5,:).*theta_s - s(3,:).*s(6,:).*theta_s;
                s(4,:).*cos(theta)      + s(3,:).*s(5,:).*theta_s - s(2,:).*s(6,:).*theta_s;
               -s(3,:).*s(4,:).*theta_s + s(5,:).*cos(theta)   + s(1,:).*s(6,:).*theta_s;
                s(2,:).*s(4,:).*theta_s - s(1,:).*s(5,:).*theta_s + s(6,:).*cos(theta);];
            % h1 = cos(theta);
            % h2 = s(1,:).*theta_s;
            % h3 = s(2,:).*theta_s;
            % h4 = s(3,:).*theta_s;
            % h5 =-s(1,:).*s(4,:).*theta_s - s(2,:).*s(5,:).*theta_s - s(3,:).*s(6,:).*theta_s;
            % h6 = s(4,:).*cos(theta)      + s(3,:).*s(5,:).*theta_s - s(2,:).*s(6,:).*theta_s;
            % h7 =-s(3,:).*s(4,:).*theta_s + s(5,:).*cos(theta)   + s(1,:).*s(6,:).*theta_s;
            % h8 = s(2,:).*s(4,:).*theta_s - s(1,:).*s(5,:).*theta_s + s(6,:).*cos(theta);
            % h = [h1;h2;h3;h4;h5;h6;h7;h8];
        end

        function x = UDQ_Log(obj,h)
            %UDQ_Log 计算h的对数映射，得到旋量形式
            % Input: h = matrix<8,n>
            %%% h = r + e*pr/2
            %%% r -> rotation -> r = cos(phi/2) + sin(phi/2)n_vec
            %%% p -> translation -> p = 0 + p_vec
            % Output: x = log(h) = phi*n_vec/2 + e*p_vec/2 -> x -> matrix<6,1>
            %%% x = phi*n
            index = find(h(1,:)==1);
            
            x1 = h(2,:).*acos(h(1,:))./sqrt(1-h(1,:).^2);
            x2 = x1.*h(3,:)./h(2,:);
            x3 = x1.*h(4,:)./h(2,:);
            x1(1,index) = 0;
            x2(1,index) = 0;
            x3(1,index) = 0;
            x4 = h(1,:).*h(6,:)-h(2,:).*h(5,:)+h(3,:).*h(8,:)-h(4,:).*h(7,:);
            x5 = h(1,:).*h(7,:) - h(2,:).*h(8,:) - h(3,:).*h(5,:) + h(4,:).*h(6,:);
            x6 = h(1,:).*h(8,:) + h(2,:).*h(7,:) - h(3,:).*h(6,:) - h(4,:).*h(5,:);
            x = [x1;x2;x3;x4;x5;x6];
        end

        function r = Q_Times(obj,varargin)
            %Q_Times 计算r1和r2的四元数乘积
            %%%Input1:  r1, r2, ... ,rn -> cell<1,n> -> matrix<4,m>
            %%%Input2:  [r1,r2, ... , rn] -> matrix<4,m>
            %%%Output:  r -> matrix<4,n>
            n = length(varargin);

            if n == 1   % 输入量为1个时，单独计算列向量相乘
                rn = varargin{1};
                m = size(rn,2); % 连乘个数
                if m<2
                    error('lack of input')
                else
                    r1 = rn(:,1);
                    for i =2:m 
                        r2 = rn(:,i);
                        % Re1 = r1(1);            Re2 = r2(1);
                        % Im1 = r1(2:4);          Im2 = r2(2:4);
                        % Re = Re1*Re2 - Im1'*Im2;
                        % Im1_hat = [0       -Im1(3)  Im1(2);
                        %            Im1(3)   0      -Im1(1);
                        %           -Im1(2)   Im1(1)  0      ];
                        % Im = Re1*Im2 + Re2*Im1 + Im1_hat*Im2;
                        % r = [Re;Im];
                        r = [r1(1)*r2(1)-r1(2)*r2(2)-r1(3)*r2(3)-r1(4)*r2(4);
                             r1(1)*r2(2)+r1(2)*r2(1)+r1(3)*r2(4)-r1(4)*r2(3);
                             r1(1)*r2(3)-r1(2)*r2(4)+r1(3)*r2(1)+r1(4)*r2(2);
                             r1(1)*r2(4)+r1(2)*r2(3)-r1(3)*r2(2)+r1(4)*r2(1);];
                        r1 = r;
                    end
                end
            else    % 输入量为多个时，分别计算对应列相乘，即n成为连乘个数，m是输出个数
                rn1 = varargin{1};
                m = size(rn1,2);
                r = zeros(4,m);
                for i =1:m  % 连乘个数
                    r1 = rn1(:,i);
                    for j =2:n % 
                        rn2 = varargin{j};
                        r2 = rn2(:,i);
                        % Re1 = r1(1);            Re2 = r2(1);
                        % Im1 = r1(2:4);          Im2 = r2(2:4);
                        % Re = Re1*Re2 - Im1'*Im2;
                        % Im1_hat = [0       -Im1(3)  Im1(2);
                        %            Im1(3)   0      -Im1(1);
                        %           -Im1(2)   Im1(1)  0      ];
                        % Im = Re1*Im2 + Re2*Im1 + Im1_hat*Im2;
                        % r(:,i) = [Re;Im];
                        r(:,i) = [r1(1)*r2(1)-r1(2)*r2(2)-r1(3)*r2(3)-r1(4)*r2(4);
                                  r1(1)*r2(2)+r1(2)*r2(1)+r1(3)*r2(4)-r1(4)*r2(3);
                                  r1(1)*r2(3)-r1(2)*r2(4)+r1(3)*r2(1)+r1(4)*r2(2);
                                  r1(1)*r2(4)+r1(2)*r2(3)-r1(3)*r2(2)+r1(4)*r2(1);];
                        r1 = r(:,i);
                    end
                end
            end
        end

        function h = UDQ_Times(obj,varargin)
            %UDQ_Times 计算h1和h2的对偶四元数乘积
            %%%Input1:  h1, h2, ... ,hn -> cell<1,n> -> matrix<8,m>
            %%%Input2:  [h1,h2, ... , hn] -> matrix<8,m>
            %%%Output:  h -> matrix<8,n>

            n = length(varargin);
            if n == 1
                hn = varargin{1};
                m = size(hn,2);
                if m <2
                    % disp('lack of input')
                    h = hn;
                else
                    h1 = hn(:,1);
                    for i = 2:m
                        h2 = hn(:,i);
                        % P1 = h1(1:4);     P2 = h2(1:4);
                        % D1 = h1(5:8);     D2 = h2(5:8);
                        % P = obj.Q_Times([P1,P2]);
                        % D = obj.Q_Times([P1,D2]) + obj.Q_Times([D1,P2]);
                        % h = [P;D];
                        
                        % P = [h1(1)*h2(1)-h1(2)*h2(2)-h1(3)*h2(3)-h1(4)*h2(4);
                        %      h1(1)*h2(2)+h1(2)*h2(1)+h1(3)*h2(4)-h1(4)*h2(3);
                        %      h1(1)*h2(3)-h1(2)*h2(4)+h1(3)*h2(1)+h1(4)*h2(2);
                        %      h1(1)*h2(4)+h1(2)*h2(3)-h1(3)*h2(2)+h1(4)*h2(1);];
                        % 
                        % D1 = [h1(1)*h2(1+4)-h1(2)*h2(2+4)-h1(3)*h2(3+4)-h1(4)*h2(4+4);
                        %       h1(1)*h2(2+4)+h1(2)*h2(1+4)+h1(3)*h2(4+4)-h1(4)*h2(3+4);
                        %       h1(1)*h2(3+4)-h1(2)*h2(4+4)+h1(3)*h2(1+4)+h1(4)*h2(2+4);
                        %       h1(1)*h2(4+4)+h1(2)*h2(3+4)-h1(3)*h2(2+4)+h1(4)*h2(1+4);];
                        % D2 = [h1(1+4)*h2(1)-h1(2+4)*h2(2)-h1(3+4)*h2(3)-h1(4+4)*h2(4);
                        %       h1(1+4)*h2(2)+h1(2+4)*h2(1)+h1(3+4)*h2(4)-h1(4+4)*h2(3);
                        %       h1(1+4)*h2(3)-h1(2+4)*h2(4)+h1(3+4)*h2(1)+h1(4+4)*h2(2);
                        %       h1(1+4)*h2(4)+h1(2+4)*h2(3)-h1(3+4)*h2(2)+h1(4+4)*h2(1);];
                        % 
                        % h = [P;D1+D2];
                        
                        h = [h1(1)*h2(1)-h1(2)*h2(2)-h1(3)*h2(3)-h1(4)*h2(4);
                             h1(1)*h2(2)+h1(2)*h2(1)+h1(3)*h2(4)-h1(4)*h2(3);
                             h1(1)*h2(3)-h1(2)*h2(4)+h1(3)*h2(1)+h1(4)*h2(2);
                             h1(1)*h2(4)+h1(2)*h2(3)-h1(3)*h2(2)+h1(4)*h2(1);

                             h1(1)*h2(1+4)-h1(2)*h2(2+4)-h1(3)*h2(3+4)-h1(4)*h2(4+4) + h1(1+4)*h2(1)-h1(2+4)*h2(2)-h1(3+4)*h2(3)-h1(4+4)*h2(4);
                             h1(1)*h2(2+4)+h1(2)*h2(1+4)+h1(3)*h2(4+4)-h1(4)*h2(3+4) + h1(1+4)*h2(2)+h1(2+4)*h2(1)+h1(3+4)*h2(4)-h1(4+4)*h2(3);
                             h1(1)*h2(3+4)-h1(2)*h2(4+4)+h1(3)*h2(1+4)+h1(4)*h2(2+4) + h1(1+4)*h2(3)-h1(2+4)*h2(4)+h1(3+4)*h2(1)+h1(4+4)*h2(2);
                             h1(1)*h2(4+4)+h1(2)*h2(3+4)-h1(3)*h2(2+4)+h1(4)*h2(1+4) + h1(1+4)*h2(4)+h1(2+4)*h2(3)-h1(3+4)*h2(2)+h1(4+4)*h2(1);];

                        h1 = h;
                    end
                end

            else    % n个h相乘
                hn1 = varargin{1};
                m = size(hn1,2); % m个列向量
                h = zeros(8,m);

                for i =2:n
                    hn2 = varargin{i};
                    % for j = 1:m
                    % h2 = hn2;
                    % P1 = hn1(1:4,:);     P2 = hn2(1:4,:);
                    % D1 = hn1(5:8,:);     D2 = hn2(5:8,:);
                    % P = obj.Q_Times(P1,P2);
                    % D = obj.Q_Times(P1,D2) + obj.Q_Times(D1,P2);
                    % h = [P;D];

                    % P = [hn1(1,:).*hn2(1,:)-hn1(2,:).*hn2(2,:)-hn1(3,:).*hn2(3,:)-hn1(4,:).*hn2(4,:);
                    %      hn1(1,:).*hn2(2,:)+hn1(2,:).*hn2(1,:)+hn1(3,:).*hn2(4,:)-hn1(4,:).*hn2(3,:);
                    %      hn1(1,:).*hn2(3,:)-hn1(2,:).*hn2(4,:)+hn1(3,:).*hn2(1,:)+hn1(4,:).*hn2(2,:);
                    %      hn1(1,:).*hn2(4,:)+hn1(2,:).*hn2(3,:)-hn1(3,:).*hn2(2,:)+hn1(4,:).*hn2(1,:);];
                    % 
                    % D1 = [hn1(1,:).*hn2(1+4,:)-hn1(2,:).*hn2(2+4,:)-hn1(3,:).*hn2(3+4,:)-hn1(4,:).*hn2(4+4,:);
                    %       hn1(1,:).*hn2(2+4,:)+hn1(2,:).*hn2(1+4,:)+hn1(3,:).*hn2(4+4,:)-hn1(4,:).*hn2(3+4,:);
                    %       hn1(1,:).*hn2(3+4,:)-hn1(2,:).*hn2(4+4,:)+hn1(3,:).*hn2(1+4,:)+hn1(4,:).*hn2(2+4,:);
                    %       hn1(1,:).*hn2(4+4,:)+hn1(2,:).*hn2(3+4,:)-hn1(3,:).*hn2(2+4,:)+hn1(4,:).*hn2(1+4,:);];
                    % 
                    % D2 = [hn1(1+4,:).*hn2(1,:)-hn1(2+4,:).*hn2(2,:)-hn1(3+4,:).*hn2(3,:)-hn1(4+4,:).*hn2(4,:);
                    %       hn1(1+4,:).*hn2(2,:)+hn1(2+4,:).*hn2(1,:)+hn1(3+4,:).*hn2(4,:)-hn1(4+4,:).*hn2(3,:);
                    %       hn1(1+4,:).*hn2(3,:)-hn1(2+4,:).*hn2(4,:)+hn1(3+4,:).*hn2(1,:)+hn1(4+4,:).*hn2(2,:);
                    %       hn1(1+4,:).*hn2(4,:)+hn1(2+4,:).*hn2(3,:)-hn1(3+4,:).*hn2(2,:)+hn1(4+4,:).*hn2(1,:);];
                    % 
                    % h = [P;D1+D2];

                    h = [hn1(1,:).*hn2(1,:)-hn1(2,:).*hn2(2,:)-hn1(3,:).*hn2(3,:)-hn1(4,:).*hn2(4,:);
                         hn1(1,:).*hn2(2,:)+hn1(2,:).*hn2(1,:)+hn1(3,:).*hn2(4,:)-hn1(4,:).*hn2(3,:);
                         hn1(1,:).*hn2(3,:)-hn1(2,:).*hn2(4,:)+hn1(3,:).*hn2(1,:)+hn1(4,:).*hn2(2,:);
                         hn1(1,:).*hn2(4,:)+hn1(2,:).*hn2(3,:)-hn1(3,:).*hn2(2,:)+hn1(4,:).*hn2(1,:);

                         hn1(1,:).*hn2(1+4,:)-hn1(2,:).*hn2(2+4,:)-hn1(3,:).*hn2(3+4,:)-hn1(4,:).*hn2(4+4,:) + hn1(1+4,:).*hn2(1,:)-hn1(2+4,:).*hn2(2,:)-hn1(3+4,:).*hn2(3,:)-hn1(4+4,:).*hn2(4,:);
                         hn1(1,:).*hn2(2+4,:)+hn1(2,:).*hn2(1+4,:)+hn1(3,:).*hn2(4+4,:)-hn1(4,:).*hn2(3+4,:) + hn1(1+4,:).*hn2(2,:)+hn1(2+4,:).*hn2(1,:)+hn1(3+4,:).*hn2(4,:)-hn1(4+4,:).*hn2(3,:);
                         hn1(1,:).*hn2(3+4,:)-hn1(2,:).*hn2(4+4,:)+hn1(3,:).*hn2(1+4,:)+hn1(4,:).*hn2(2+4,:) + hn1(1+4,:).*hn2(3,:)-hn1(2+4,:).*hn2(4,:)+hn1(3+4,:).*hn2(1,:)+hn1(4+4,:).*hn2(2,:);
                         hn1(1,:).*hn2(4+4,:)+hn1(2,:).*hn2(3+4,:)-hn1(3,:).*hn2(2+4,:)+hn1(4,:).*hn2(1+4,:) + hn1(1+4,:).*hn2(4,:)+hn1(2+4,:).*hn2(3,:)-hn1(3+4,:).*hn2(2,:)+hn1(4+4,:).*hn2(1,:)];
                    hn1 = h;
                    % end
                end
            end
        end

        function h_s = UDQ_Conj(obj,h)
            %UDQ_CONJ(h)计算h的共轭
            %%%Input:   h -> matrix<8,n>
            %%%Output:  h_s -> matrix<8,n>
            n = size(h,2);
            ReH = zeros(8,n);
            ReH([1,5],:) = h([1,5],:);
            ImH = h; ImH([1,5],:) = 0;
            h_s = ReH - ImH;
        end

        function h_inv = UDQ_Inv(obj,h)
            %UDQ_INV(h) 计算h的逆
            %%%Input:   h -> matrix<8,n>
            %%%Output:  h_inv -> matrix<8,n>
            n = size(h,2);
            h_s = obj.UDQ_Conj(h);
            check = obj.UDQ_Times(h,h_s);
            check (abs(check)<1e-5) = 0;
            check (abs(check-1)<1e-5)=1;
            judge = isequal(check, repmat([1;0;0;0;0;0;0;0],1,n));
            if ~judge
                error('something wrong')
            else
                h_inv = h_s ;
            end
        end

        function Ad_h = UDQ_Ad(obj,h)
            %UDQ_AD(h) 计算UDQ h:SE(3) 的伴随表示
            %%%Input:   h -> matrix<8,1> -> SE(3) -> spin(3) HDP R3
            %%%Output:  Ad_h -> matrix<6,6> -> se(3) -> 
            r11 = h(1)^2 + h(2)^2  - h(3)^2 - h(4)^2;
            r12 = 2*(h(2)*h(3) - h(1)*h(4));
            r13 = 2*(h(1)*h(3) + h(2)*h(4));
            r21 = 2*(h(2)*h(3) + h(1)*h(4));
            r22 = h(1)^2 - h(2)^2  + h(3)^2 - h(4)^2;
            r23 = 2*(h(3)*h(4) - h(1)*h(2));
            r31 = 2*(h(2)*h(4) - h(1)*h(3));
            r32 = 2*(h(1)*h(2) + h(3)*h(4));
            r33 = h(1)^2 - h(2)^2  - h(3)^2 + h(4)^2;
            R = [r11,r12,r13;
                r21,r22,r23;
                r31,r32,r33;];

            m11 = 2*(h(1)*h(5) + h(2)*h(6) - h(3)*h(7) - h(4)*h(8));
            m12 = 2*(h(2)*h(7) + h(3)*h(6) - h(4)*h(5) - h(1)*h(8));
            m13 = 2*(h(1)*h(7) + h(2)*h(8) + h(3)*h(5) + h(4)*h(6));
            m21 = 2*(h(1)*h(8) + h(2)*h(7) + h(3)*h(6) + h(4)*h(5));
            m22 = 2*(h(1)*h(5) - h(2)*h(6) + h(3)*h(7) - h(4)*h(8));
            m23 = 2*(h(3)*h(8) + h(4)*h(7) - h(1)*h(6) - h(2)*h(5));
            m31 = 2*(h(2)*h(8) - h(3)*h(5) + h(4)*h(6) - h(1)*h(7));
            m32 = 2*(h(1)*h(6) + h(2)*h(5) + h(3)*h(8) + h(4)*h(7));
            m33 = 2*(h(1)*h(5) - h(2)*h(6) - h(3)*h(7) + h(4)*h(8));
            M = [m11,m12,m13;
                m21,m22,m23;
                m31,m32,m33;];
            Ad_h = zeros(6);
            Ad_h(:,:) = [R,zeros(3);M,R];
        end

        function Z = Twist_Ad(obj,t)
            %TWIST_AD(t) 计算旋量t：se3 的伴随表示
            %%%Input:   t -> matrix<6,1>
            %%%Output:  Z -> matrix<6,6>
            Z = [0,    -t(3),   t(2),   0,      0,      0;
                 t(3),  0,     -t(1),   0,      0,      0;
                -t(2),  t(1),   0,      0,      0,      0;
                 0,    -t(6),   t(5),   0,     -t(3),   t(2);
                 t(6),  0,     -t(4),   t(3),   0,     -t(1);
                -t(5),  t(4),   0,     -t(2),   t(1),   0;];
            Z = 2*Z;
        end

        function At = Twist_Jacobian(obj,t)
            %TWIST_JACOBIAN(t) 计算旋量t的雅可比矩阵
            %%%Input:   t -> matrix<6,1>
            %%%Output:  At -> matrix<6,6>
            Z = obj.Twist_Ad(t);
            phi = sqrt(t(1)^2 + t(2)^2 + t(3)^2)*2;
            if phi == 0 
                At = eye(6) + Z/2;
            else
                At = eye(6) + (4-phi*sin(phi)-4*cos(phi))    / (2*phi^2) * Z + ...
                              (4*phi-5*sin(phi)+phi*cos(phi))/ (2*phi^3) * Z^2 + ...
                              (2-phi*sin(phi)-2*cos(phi))    / (2*phi^4) * Z^3 + ...
                              (2*phi-3*sin(phi)+phi*cos(phi))/ (2*phi^5) * Z^4 ;
            end
        end

        function r = Q_Log2r(obj,R)
            %Q_LOG2R(R) 计算旋转矩阵R的四元数表示
            %%%Input:   R -> matirx<3,3>
            %%%Output:  r -> matrix<4,1>
            traceR=sum(diag(R));
            tmp=(traceR-1)/2; % tr(R) = 1+2cos(theta)
            if tmp>1
                tmp=1;
            elseif tmp<-1
                tmp=-1;
            end

            theta=acos(tmp);    % 返回值为[0,pi]
            if theta==0            %theta=0 轴w任意
                w=[0.0;0.0;0.0];
            elseif theta-pi==0     %theta=pi
                w=null(eye(3)-R); %这是因为(I-R)w=0 等效于w=Rw，旋转矩阵作用于旋转轴向量不会起作用，相当于一个不动点
                w=w(:,1);%选取第一列当作旋转轴
                nw=abs(w);
                if nw(1)>0.1
                    w=w(1)/nw(1)*w;
                elseif nw(2)>0.1
                    w=w(2)/nw(2)*w;
                else
                    w=w(3)/nw(3)*w;
                end
                w=w*pi;
            else
                w=[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)]/(2*sin(theta))*theta;
            end

            phi = norm(w);
            r = [cos(phi/2);sin(phi/2)*w/phi];
            if phi < 1e-9
                r = [cos(phi/2);w/2];
            end
        end

        function [h_a,h_n,Q_n0] = Cail_Pose2UDQ(obj,g_a,g_n,q_n)
            %CAIL_POSE2UDQ(g_a, g_n, q_n) 将标定所需数据从 HTM 转化为 UDQ
            %%%Input:   g_a -> matrix<4,4,n>, g_n -> matrix<4,4,n>
            %%%         q_n -> matrix<6,n> -> 各支链主动关节
            %%%Output:  h_a -> matrix<8,n>, h_n -> matrix<8,n>
            %%%         Q_n0-> matrix<6,6,n> -> 各支链各关节
            n = size(g_a,3);
            p_n = zeros(4,n); p_n(2:4,:) = g_n(1:3,4,:);
            p_a = zeros(4,n); p_a(2:4,:) = g_a(1:3,4,:);
            r_n = zeros(4,n);
            r_a = zeros(4,n);
            Q_n0 = zeros(6,6,n);
            for i =1:n
                r_n(:,i) = obj.Q_Log2r(g_n(1:3,1:3,i));
                r_a(:,i) = obj.Q_Log2r(g_a(1:3,1:3,i));
                Q_n0(3,:,i) = q_n(:,i);
            end
            h_n = [r_n;obj.Q_Times(p_n,r_n)/2];
            h_a = [r_a;obj.Q_Times(p_a,r_a)/2];
        end
        
        function [h] = Pose2UDQ(obj,g)
            %POSE2UDQ(g) 将位姿g（HTM）转化为h（UDQ）
            %%%Input:   g -> matrix<4,4,n>
            %%%Output:  h -> matrix<8,n>
            n = size(g,3);
            p = zeros(4,n); p(2:4,:) = g(1:3,4,:);
            r = zeros(4,n);
            for i =1:n
                r(:,i) = obj.Q_Log2r(g(1:3,1:3,i));
            end
            h = [r;obj.Q_Times(p,r)/2];
        end

        function [g] = UDQ2Pose(obj,h)
            %UDQ2POSE(h) 将UDQ形式位姿转化为HTM位姿
            %%%Input:   h -> matrix<8,n>
            %%%Output:  g -> matrix<4,4,n>
            n = size(h,2);
            g = zeros(4,4,n);
            g(4,4,:) = ones(1,n);
            [r,~] = obj.PD(h);
            % check = obj.Q_Times(r,[r(1);-r(2:4)])
            for i = 1:n
                % theta = asin(r(1,i))*2;
                R = [1-2*r(3,i)^2-2*r(4,i)^2, 2*r(2,i)*r(3,i)-2*r(4,i)*r(1,i), 2*r(2,i)*r(4,i)+2*r(3,i)*r(1,i);
                     2*r(2,i)*r(3,i)+2*r(4,i)*r(1,i), 1-2*r(2,i)^2-2*r(4,i)^2, 2*r(3,i)*r(4,i)-2*r(2,i)*r(1,i);
                     2*r(2,i)*r(4,i)-2*r(3,i)*r(1,i), 2*r(3,i)*r(4,i)+2*r(2,i)*r(1,i), 1-2*r(2,i)^2-2*r(3,i)^2];
                g(1:3,1:3,i) = R;
            end
            h_s = obj.UDQ_Conj(h(:,i));
            p0 = obj.Q_Times(2*h(5:8,:),h_s(1:4,:));
            g(1:3,4,:) = p0(2:4,:);
            
        end

        function [pe] = UDQ_err(obj,ha,hn)
            num = size(ha,2);
            if num ~= size(hn,2)
                error('输入量不匹配')
            end
            pe_total = zeros(4,num);
            % we = zeros(3,num);
            for i = 1:num
                % r+ e 0.5pr
                % ra = ha(1:4,:);
                ha_s = obj.UDQ_Conj(ha);
                hn_s = obj.UDQ_Conj(hn);
                pa = 2*obj.Q_Times(ha(5:8,:),ha_s(1:4,:));
                pn = 2*obj.Q_Times(hn(5:8,:),hn_s(1:4,:));
                pe_total = pa-pn;
            end
            pe = vecnorm(pe_total(2:4,:));
        end

        function X = UDQ_LeastSquare(obj,A,B,X0)
            % solve AX = B least square problem
            % A in Q{m,n}, X in Q{n,p}, B in Q{m,p}
            Z = X0;
            % problem 1. A 的结构
            % problem 2. Ast 的伪逆
            X1 = A1\B1 + (eye(n)-A1\A1)*Z;
            C = A2*(eye(n)-A1\A1);
            D = A2*A1\B1 - B2;
            G = [A1,C];
            W1 = [X2,Z];
            W1k = W1;
            tau = 0.5;
            r = rank(G);
            for num = 1:50
                W_temp = W1k;
                W1k = (tau*eye(r) + G'*G)\(tau*W1k-G1'*D);
                if (norm(W1k - W_temp)<1e-9)
                    break
                end
            end
            Wk = V'*[W1k;W2k];
            X2k = Wk(1:n,:);
            Zk = Wk(n+1:2*n,:);
            X1k = A1\B1+(eye(n)-A1\A1)*Zk;
            X = [X1k;X2k];

        end

        % 
        function [P,D,W] = Q_SVD4(AA)
            % This is a real structure-preserving algorithm for QSVD
            % A is in Q{m,n}
            % A = V*S*U^H
            B4 = AA;
            P = zeros(4*m, m);
            P(1 : m, :) = eye(m);
            W = zeros(4*n, n);
            W (1 : n, :) = eye(n);
            Y = W';
            for s = 1 : n - 1 
                if norm([B4(s:m,s); B4((m+s):(2*m),s); B4((2*m+s):(3*m),s); B4((3*m+s):(4*m),s)]) ~= 0 
                    if norm([B4(s,s); B4(s+m,s); B4(s+2*m,s); B4(s+3*m,s)]) ~= 0
                        G4 = JRSGivens(B4(s,s), B4(s+m,s), B4(s+2*m,s), B4(s+3*m, s));
                        t = s;
                        B4([t,t+m,t+2*m,t+3*m], s:n) = G4'*B4([t,t+m,t+2*m,t+3*m], s:n);
                        P([t, t+m, t+2*m, t+3*m], 1:m) = G4'*P([t,t+m,t+2*m,t+3*m], 1:m); 
                    end 
                    [u4R, beta4] = Householder(B4(s:m,s),B4((m+s):(2*m),s),B4((2*m+s):(3*m),s),B4((3*m+s):(4*m),s),m-s+1);
                    BB = B4([s:m,s+m:2*m,s+2*m:3*m,s+3*m:4*m],s:n);
                    BB = BB-(beta4*u4R)*(u4R'*BB);
                    B4([s:m, s+m:2*m, s+2*m:3*m, s+3*m:4*m], s:n) = BB;
                    PP = P([s:m, s+m:2*m, s+2*m:3*m,s+3*m:4*m],1:m);
                    PP = PP-(beta4*u4R)*(u4R'*PP);
                    P([s:m,s+m:2*m,s+2*m:3*m,s+3*m:4*m],1:m)=PP;
                end 
                if s <= n-2
                    if norm([B4(s, s + 1 : n)';-B4(s+m,s+1:n)';-B4(s+2*m,s+1:n)';-B4(s+3*m,s+1:n)']) ~= 0
                        Z(s:m, [s+1:n, s+1+n:2*n, s+1+2*n:3*n, s+1+3*n:4*n]) = [B4(s:m,s+1:n),-B4(s+m:2*m, s+1:n),-B4(s+2*m:3*m, s+1:n),-B4(s+3*m:4*m, s+1:n)];
                        if norm([Z(s,s+1);Z(s,s+1+n);Z(s,s+1+2*n);Z(s,s+1+3*n)]) ~= 0
                            G4 = JRSGivens(Z(s,s+1),Z(s,s+1+n), Z(s,s+1+2*n),Z(s,s+1+3*n));
                            t=s+1;
                            Z(s:m, [t,t+n, t+2*n, t+3*n]) = Z(s:m,[t,t+n,t+2*n, t+3*n])*G4;
                            Y(1:n, [t,t+n, t+2*n, t+3*n]) = Y (1:n, [t,t+n,t+2*n,t+3*n])*G4;
                        end 
                        [u4R, beta4] = Householder(Z(s,s+1:n)',Z(s,s+1+n:2*n)',Z(s,s+1+2*n:3*n)', Z(s,s+1+3*n:4*n)', n-s);
                        ZZ = Z(s:m, [s+1:n,s+1+n:2*n,s+1+2*n:3*n,s+1+3*n:4*n]);
                        ZZ = ZZ - ZZ*u4R*(beta4*u4R');
                        ns = n-s;
                        B4([s:m, s+m:2*m, s+2*m:3*m, s+3*m:4*m],s+1:n) = [ZZ(:,1:ns);-ZZ(:,ns+1:2*ns);-ZZ(:,2*ns+1:3*ns);-ZZ(:,3*ns+1:4*ns)];
                        YY = Y(1:n,[s+1:n,s+1+n:2*n, s+1+2*n:3*n, s+1+3*n:4*n]);
                        YY = YY*YY*u4R*(beta4*u4R');
                        Y(1:n,[s+1:n,s+1+n:2*n,s+1+2*n:3*n,s+1+3*n:4*n]) = YY ;
                    end
                elseif s == n - 1
                        if norm([B4(s:m,s+1:n);-B4(s+m:2*m,s+1:n);-B4(s+2*m:3*m,s+1:n);-B4(s+3*m:4*m,s+1:n)]) ~= 0 
                            Z(s:m,[s+1:n,s+1+n:2*n,s+1+2*n:3*n,s+1+3*n:4*n]) = [B4(s:m,s+1:n), -B4(s+m:2*m,s+1:n),-B4(s+2*m:3*m,s+1:n),-B4(s+3*m:4*m,s+1:n)];
                            G4 = JRSGivens(Z(s,s+1), Z(s,s+1+n), Z(s,s+1+2*n), Z(s,s+1+3*n));
                            t = s+1;
                            Z(s:m,[t,t+n,t+2*n,t+3*n]) = Z(s:m, [t,t+n, t+2*n, t+3*n])*G4;
                            Y (1:n, [t,t+n,t+2*n, t+3*n]) = Y (1:n, [t,t+n, t+2*n, t+3*n])*G4;
                            B4([s:m, s+m:2*m, s+2*m:3*m, s+3*m:4*m], s+1:n) = [Z(s:m, s+1:n);-Z(s:m,s+1+n:2*n);-Z(s:m,s+1+2*n:3*n);-Z(s:m,s+1+3*n:4*n)];
                        end
                end
            end
            s = n;
            if norm([B4(s:m,s);B4((m+s):(2*m),s); B4((2*m+s):(3*m),s); B4((3*m+s):(4*m),s)]) ~= 0 
                if norm([B4(s,s); B4(s+m,s);B4(s+2*m,s);B4(s+3*m,s)]) ~= 0
                    G4 = JRSGivens(B4(s,s),B4(s+m,s), B4(s+2*m,s), B4(s+3*m,s));
                    t = s;
                    B4([t,t+m,t+2*m,t+3*m],s:n) = G4'*B4([t,t+m,t+2*m,t+3*m],s:n); 
                    P([t,t+m, t+2*m, t+3*m], 1:m) = G4'*P([t,t+m,t+2*m,t+3*m],1:m);
                end
                if m > n
                    [u4R, beta4] = Householder(B4(s:m,s),B4((m+s):(2*m),s),B4((2*m+s):(3*m),s),B4((3*m+s):(4*m),s),m-s+1);
                    BB4 = B4([s:m,s+m:2*m,s+2*m:3*m,s+3*m:4*m],s:n);
                    BB4 = BB4 -(beta4*u4R)*(u4R'*BB4);
                    B4([s:m,s+m:2*m,s+2*m:3*m, s+3*m:4*m], s:n) = BB4;
                    PP4 = P([s:m,s+m:2*m,s+2*m:3*m,s+3*m:4*m],1:m);
                    PP4 = PP4 - (beta4-u4R)*(u4R'*PP4);
                    P([s:m, s+m:2*m, s+2*m:3*m, s+3*m:4*m], 1:m) = PP4;
                end
            end
            [U, D,V] = svd(B4(1:m,1:n));
            P = [U'*P(1:m,:);U'*P(m+1:2*m, :);U'*P(2*m+1:3*m,:);U'*P(3*m+1:4*m,:)];
            P = [P(1:m,:)';-P(m+1:2*m,:)';-P(2*m+1:3*m,:)';-P(3*m+1:4*m, :)'];
            Y = [Y(:,1:n)*V, Y(:,n+1:2*n)*V, Y(:,2*n+1:3*n)*V, Y(:,3*n+1:4*n)*V];
            W = [Y(:,1:n);-Y(:,n+1:2*n);-Y(:,2*n+1:3*n);-Y(:,3*n+1:4*n)];
        end

        function G2 = JRSGivens(g1,g2,g3,g4)
            % g = g0 + g1i + g2j + g3k
            if [g2,g3,g4] == 0
                G2 = eye(4);
            else
                G2 = Realp(g1,g2,g3,g4)/norm([g1,g2,g3,g4]);
            end
        end

        function GR = Realp(g1,g2,g3,g4)
            GR = [g1,-g2,-g3,-g4;
                  g2, g1,-g4, g3;
                  g3, g4, g1,-g2;
                  g4,-g3, g2, g1];
        end

        function [u,b] = Householder0(obj,x,n)
            % Real Householder Transformation
            % u is Householder vector with H(u) = I - 2*u*u'
            % beta is the scalar with H = I - beta*x*x'
            % x is in R^n
            u = x;
            a1 = norm(x);
            % u = y-x, Hx = y
            if u(1)>=0
                u(1) = u(1)+a1;
                b = 1/(u(1)*a1);
            else
                u(1)=u(1)-a1;
                b=1/(u(1)*a1);
            end
        end

        function  [u,b1] = Householder1(obj,x1,x2,x3,x4,n)
            % Quaternion Householder Based Transformation
            % H1
            % x = x1 + x2i + x3j + x4k
            % n 表示输入四元数个数
            % u 是四元数 u1 的实矩阵表示
            
            u1(1:n,1:4) = [x1,x2,x3,x4];
            aa = norm([x1;x2;x3;x4]);
            xx = norm([x1(1),x2(1),x3(1),x4(1)]);
            % 计算映射后的四元数 a1
            if xx == 0
                a1 = aa*[1,0,0,0];
            else
                a1 = -(aa/xx)*([x1(1),x2(1),x3(1),x4(1)]);
            end

            u1(1,1:4) = u1(1,1:4)-a1;
            b1 = 1/(aa*(aa+xx));
            u = Realp(u1(:,1),u1(:,2),u1(:,3),u1(:,4));
        end
        
        function [u,b] = Householder4(x1,x2,x3,x4,n)
            H1 = Householder1(x1,x2,x3,x4,n);
            
        end
        
        function xt = ScLerp(obj,xa,xb,t)
            n = length(t);
            xt = zeros(8,n);
            xa_star = obj.UDQ_Conj(xa);
            for i = 1:n
                xexp = obj.UDQ_Exp(t(i)*obj.UDQ_Log(obj.UDQ_Times(xa_star,xb)));
                xt(:,i) = obj.UDQ_Times(xa,xexp);
            end
        end

        function q = UDQ_rand(obj,n)
            q = zeros(8,n);
            w = rand(3,n)*10-5;
            hat = @(x)[0   ,-x(3), x(2);
                       x(3), 0   ,-x(1);
                      -x(2), x(1), 0];
            for i = 1:n
                w_hat = hat(w(:,i));
                R = expm(w_hat);
                p = rand(3,1)*20-10;
                g = [R,p;0,0,0,1];
                q(:,i) = obj.Pose2UDQ(g);
            end
        end
    end
end