function [p,R] = fkine(a,b,Len,a_drive,n)
% forward kinematics
% 对偶四元数方法

% [w,theta] = Quaternion(R);

% Representation of pose of mobile platform


conj = @(x)[-x(1:3);x(4)];
skew = @(x)[0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
% E_1 = p + Rb =  d + l_i;
% E_2(theta)*xi = lambda + xi*ai + ci*xi
% norm(E_2)^2 = **** = norm(l_i)^2 = x'Qx/2 - C_i

Q = zeros(8,8,8);           
Q(1:4, 1:4, 7)=2*eye(4); % xi 模长为1 
Q(1:4, 5:8, 8)=eye(4);
Q(5:8, 1:4, 8)=eye(4);   % xi 和 lambda 正交 


C = [Len.^2,1,0]';

a(2,:) = a_drive;

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
xi_0 = [0 0 0 1]';           % 四元数-旋转相关
lambda_0 = [0 -1 0 0]';      % 四元数-平移相关
X = [xi_0;lambda_0];

while norm(flag)>1e-10
    J = [Q(:,:,1)*X, Q(:,:,2)*X, Q(:,:,3)*X, Q(:,:,4)*X, Q(:,:,5)*X, Q(:,:,6)*X, Q(:,:,7)*X, Q(:,:,8)*X]';
    DET_J = det(J);
    % % delta_X = (J'*J)\(J'*C);
    delta_X = pinv(J)*C;
    % delta_X = J\C;
    X1 = 0.5*X + delta_X;
    flag =  X1-X;
    X=X1;   
    
    % norm(flag)
end

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
p = quaternion_multi(lambda,conj(xi))*n;
% R = rotation_matrix(w,theta);
R = skew(w)*sin(theta) + skew(w)*skew(w)*(1-cos(theta)) + eye(3);