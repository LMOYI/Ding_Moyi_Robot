% 适用于3维列向量叉乘
function qx = CrossMatrix(q)
qx = zeros(3);
qx(1,2) = -q(3);
qx(2,1) =  q(3);
qx(1,3) =  q(2);
qx(3,1) = -q(2);
qx(2,3) = -q(1);
qx(3,2) =  q(1);
end