function [a,b]=Initial_Config(A_d, A_L, B_d, B_L,k,n,varargin)
%INTIAL_CONFIG [a,b]=Initial_Config(A_d, A_L, B_d, B_L,k,n)
%   Input : n 控制输出单位，其他尺寸输入均按照mm来。
%   Output: a,b 单位按照n来
if nargin < 6
    n = 1;
end
% 右上角起，顺时针排序
if k == -1
    ra = varargin{1}/2;
    rb = varargin{2}/2;

    theta_a1 = acos(A_L/2/ra);
    theta_a2 = asin(A_d/2/ra);
    a1 = theta_a1;
    a2 = theta_a1 + 2*theta_a2;
    a3 = theta_a2;

    theta_b1 = acos(A_L/2/rb);
    theta_b2 = asin(A_d/2/rb);
    b1 = theta_b1;
    b2 = theta_b1 - 2*theta_b2;
    b3 = theta_b2;

    a = [ra*cos(a2), 0, ra*sin(a2);
         ra*cos(a1), 0, ra*sin(a1);
         ra*sin(a3), 0,-ra*cos(a3);
        -ra*sin(a3), 0,-ra*cos(a3);
        -ra*cos(a1), 0, ra*sin(a1);
        -ra*cos(a2), 0, ra*sin(a2);]'/n;

    b = [rb*sin(b3), 0, rb*cos(b3);
         rb*cos(b2), 0,-rb*sin(b2);
         rb*cos(b1), 0,-rb*sin(b1);
        -rb*cos(b1), 0,-rb*sin(b1);
        -rb*cos(b3), 0, rb*sin(b3);
        -rb*sin(b2), 0, rb*cos(b2);]'/n;
elseif k == 0       % 方案0
    theta1 = atan(A_d/A_L);
    r1 = sqrt(A_d^2+A_L^2)/2;
    theta2 = atan(B_d/B_L);
    r2 = sqrt(B_d^2+B_L^2)/2;
    a = [r1*cos(theta1), 0, r1*sin(theta1);
         r1*cos(theta1), 0,-r1*sin(theta1);
         r1*sin(theta1), 0,-r1*cos(theta1);
        -r1*sin(theta1), 0,-r1*cos(theta1);
        -r1*cos(theta1), 0,-r1*sin(theta1);
        -r1*cos(theta1), 0, r1*sin(theta1);]'/n;
    
    b = [r2*cos(theta2), 0, r2*sin(theta2);
         r2*cos(theta2), 0,-r2*sin(theta2);
         r2*sin(theta2), 0,-r2*cos(theta2);
        -r2*sin(theta2), 0,-r2*cos(theta2);
        -r2*cos(theta2), 0,-r2*sin(theta2);
        -r2*cos(theta2), 0, r2*sin(theta2);]'/n;

elseif k == 1
    a = [A_L/2, 0, B_L/2;
         A_L/2, 0, B_L/2-A_d;
         A_d/2, 0,-A_L/2;
        -A_d/2, 0,-A_L/2;
        -A_L/2, 0, B_L/2-A_d;
        -A_L/2, 0, B_L/2]'/n;
    b = [B_L/2, 0, B_L/2;
         B_L/2, 0, B_L/2-B_d;
         B_d/2, 0,-B_L/2;
        -B_d/2, 0,-B_L/2;
        -B_L/2, 0, B_L/2-B_d;
        -B_L/2, 0, B_L/2]'/n;

elseif k==2
    a = zeros(3,6);
    a(:,1:3) = [A_L/2, 0, A_d;
                A_L/2, 0, 0;
                A_d/2, 0,-A_L/2;]'/n;
    a(:,4:6) = a(:,[3,2,1]);
    a(1,4:6) = -a(1,4:6);

    b = zeros(3,6);
    b(:,1:3) = [B_L/2, 0, B_d;
                B_L/2, 0, 0;
                B_d/2, 0, -B_L/2]'/n;
    b(:,4:6) = b(:,[3,2,1]);
    b(1,4:6) = -b(1,4:6);

else
    a = [A_d/2, 0, A_L/2;
        A_d/2, 0, 0;
        A_d/2, 0,-A_L/2;
        -A_d/2, 0,-A_L/2;
        -A_d/2, 0, 0;
        -A_d/2, 0, A_L/2;]'/n;
    b = [B_d/2, 0, B_L/2;
        B_d/2, 0, 0;
        B_d/2, 0,-B_L/2;
        -B_d/2, 0,-B_L/2;
        -B_d/2, 0, 0;
        -B_d/2, 0, B_L/2;]'/n;
end