function [Jac,Jac_inv] = iJacobian(a,b,Len,p,R) 
% S^dot = J^inv * V , V = [v;w]
% Jac 为几何雅可比
p = p(1:3); 
Jac_inv = zeros(6,6);
for i = 1:6
    g = R*b(:,i);
    rho = sqrt(Len(i)^2-(p(1)+g(1)-a(1,i))^2-(p(3)+g(3)-a(3,i))^2 ); % 连杆向量y向分量
    r =  -[p(1)+g(1)-a(1,i), -rho, p(3)+g(3)-a(3,i)]';   %
    G =  skew(g);
    Jac_inv(i,:) = [r', (G*r)']/rho;
end
Jac = inv(Jac_inv);
end