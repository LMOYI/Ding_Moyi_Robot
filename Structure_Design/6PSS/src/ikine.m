function S = ikine(a,b,Len,p,R)

p = p(1:3);
S = zeros(1,6); % P副上球铰相对位置（在y向）-反映驱动距离

for i = 1:6
B = p+R*b(:,i);
S(i) = B(2)+sqrt(Len(i)^2 - (B(1)-a(1,i))^2 - (B(3)-a(3,i))^2);
end
end