classdef mySE3
    %SE3 This file is a module for calculation with SE(3) and se(3).
    %   Author: Weijia Zhang.
    %   version 1.0
    
    properties
    %该类没有定义任何属性，所有操作都是通过方法实现的。
    end
    
    methods
        function obj = mySE3()
            %SE3 构造此类的实例
        end
        
        function S = Sw(obj,w)
            %Sw 求向量对应的反对称矩阵
            %   input:
            %   w::matrix<3,1>
            %   output:
            %   S::matrix<3,3>
            S=[0,-w(3),w(2);
                w(3),0,-w(1);
                -w(2),w(1),0];
        end
        
        function R = exp2R(obj,w)
            %exp2R 求向量对应的旋转矩阵
            %通过向量的方向表示旋转轴，向量的模长表示旋转角度，应用罗德里格斯公式计算旋转矩阵 
            %   input:
            %   w::matrix<3,1>
            %   output:
            %   R::matrix<3,3>
            theta=norm(w);
            if theta==0
                R=eye(3);
            else
                w1=w/theta;
                W=obj.Sw(w1);
                R=eye(3)+W*sin(theta)+W*W*(1-cos(theta));
            end
        end
        
        function T = exp2T(obj,s)
            %exp2T 求向量对应的齐次变换矩阵
            %输入为旋量，3.22ppt有转换公式，其中v需要注意除去theta
            %   input:
            %   s::matrix<6,1>
            %   output:
            %   T::matrix<4,4>
            w=s(1:3);
            v=s(4:6);
            R=obj.exp2R(w);
            theta=norm(w);
            if theta<1e-15
                T=[[R v]; [0 0 0 1.0]];
            else
                w=w/theta;
                v=v/theta;
                W=obj.Sw(w);
                T=[[R (eye(3)-R)*W*v+(w'*v)*w*theta]; [0 0 0 1.0]];
            end
        end
        
        function w = log2w(obj,R)
            %log2w 旋转矩阵的对数运算
            %
            %   input:
            %   R::matrix<3,3>
            %   output:
            %   w::matrix<3,1>
            traceR=sum(diag(R));
            tmp=(traceR-1)/2;
            if tmp>1
                tmp=1;
            elseif tmp<-1
                tmp=-1;
            end
            theta=acos(tmp);
            if theta==0            %theta=0
                w=[0.0;0.0;0.0];
            elseif theta-pi==0      %theta=pi
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
                w=[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)]/(2*sin(theta));
                w=w*theta;
            end
        end
        
        function s = log2s(obj,T)
            %log2s 齐次变换矩阵的对数运算  矩阵转换为旋量
            %   input:
            %   T::matrix<4,4>
            %   output:
            %   s::matrix<6,1>
            R=T(1:3,1:3);
            p=T(1:3,4);
            w=obj.log2w(R);
            theta=norm(w);
            if theta<1e-16
                v=p;
                s=[[0,0,0]'; v];
            else
                w=w/theta;
                W=obj.Sw(w);
                A=(eye(3)-R)*W+w*w'*theta;
                if rank(A)<3
                    v=(A+0.00001*eye(3))\p;
                else
                    v=A\p;
                end
                s=[w; v]*theta;
            end
        end
        
        function A = Ad(obj,T)
            %Ad 求SE3的伴随矩阵
            %   input:
            %   T::matrix<4,4>
            %   output:
            %   A::matrix<6,6>
            R=T(1:3,1:3);
            p=T(1:3,4);
            Sp=obj.Sw(p);
            A=[[R zeros(3,3)]; [Sp*R R]];
        end
        
        function A = ad(obj,s)
            %ad 求se3的伴随矩阵
            %   input:
            %   s::matrix<6,1>
            %   output:
            %   A::matrix<6,6>
            w=s(1:3);
            v=s(4:6);
            A=[[obj.Sw(w) zeros(3,3)]; [obj.Sw(v) obj.Sw(w)]];
        end
        
        function Tv = Tinv(obj,T)
            %Tinv 齐次变换矩阵求逆
            %   input:
            %   T::matrix<4,4>
            %   output:
            %   Tv::matrix<4,4>
            R=T(1:3,1:3);
            p=T(1:3,4);
            Tv=[[R' -R'*p]; [0 0 0 1.0]];
        end
        
        function R = rotx(obj,theta)
            %rotx 沿x方向的旋转矩阵
            %   input:
            %   theta::double
            %   output:
            %   R::matrix<3,3>
            w=[1.0, 0, 0]';
            R=obj.exp2R(w*theta);
        end
        
        function R = roty(obj,theta)
            %rotx 沿y方向的旋转矩阵
            %   input:
            %   theta::double
            %   output:
            %   R::matrix<3,3>
            w=[0, 1.0, 0]';
            R=obj.exp2R(w*theta);
        end
        
        function R = rotz(obj,theta)
            %rotx 沿z方向的旋转矩阵
            %   input:
            %   theta::double
            %   output:
            %   R::matrix<3,3>
            w=[0, 0, 1.0]';
            R=obj.exp2R(w*theta);
        end
        
        function p=T2pose(obj, T)
            % this funciton is for switch the transformation matrix to
            % the defaut pose form [x,y,z,rx,ry,rz]'
            % unit: mm, rad
            p=T(1:3,4);
            R=T(1:3,1:3);
            rot=rotm2eul(R,'ZYX');%ZYX,动轴，右乘
            rx=rot(3);
            %rx=rx/pi*180;
            ry=rot(2);
            %ry=ry/pi*180;
            rz=rot(1);
            %rz=rz/pi*180;
            p = [p;rx;ry;rz];
        end

        function T=gst2F(gst)
            % this function is to changethe gst of [x,y,z,rx,ry,rz] to frame of end
            % relative to the base frame
            p=gst(1:3);
            p=reshape(p,[3 1]);
            R = mySE3.rotz(gst(6))*mySE3.roty(gst(5))*mySE3.rotx(gst(4));
            T = [R,p;0,0,0,1];
        end
    end
end

