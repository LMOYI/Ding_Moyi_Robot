function [drive_force,plat_force] = force(a,b,s,p,R,F,T)
% F = [0 0 -1200]';   % 末端载荷
% T = [0 100 200]';   % 末端扭矩

L = zeros(3,6);     % 支链向量
Lx = zeros(3,3,6);  % 支链向量叉乘矩阵
% bx = zeros(3,3,6);  % 末端铰链叉乘矩阵（用于计算动平台力矩平衡）
Bx = zeros(3,3,6);  % 末端铰链叉乘矩阵（用于计算支链力矩平衡）
m = 0;              % 连杆质量 kg
M = zeros(3,6); M(3,:) = -m*10*ones(1,6); % 质量矩阵
N = zeros(3,6);     % 质量引起力矩（参考坐标系下，各支链S1为支点）
a(2,:) = s;
for i =1:6
    L(:,i) = p + R*b(:,i) - a(:,i);
    Lx(:,:,i) = CrossMatrix(L(:,i));
    % bx(:,:,i) = CrossMatrix(R*b(:,i));
    N(:,i) = Lx(:,:,i)/2*M(:,i); 
    % Bx(:,:,i) = CrossMatrix(p+R*b(:,i));
    Bx(:,:,i) = CrossMatrix(R*b(:,i));
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
% eq_A*[F_a; F_b]=eq_B
eq_B = [-reshape(M,[18,1]);
        -reshape( ...
        ...
        ...
        N,[18,1]);
        F;
        % T+CrossMatrix(p)*F
        T];
% rank(eq_A)-rank([eq_A,eq_B])
x_force = eq_A\eq_B;

drive_force = reshape(x_force(1:18),[3,6]);
plat_force  = reshape(x_force(19:end),[3,6]);
end
