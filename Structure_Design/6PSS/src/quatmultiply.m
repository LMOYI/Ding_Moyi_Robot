function qr=quatmultiply(q,r)
qr = zeros(4,1);
qr(1) = r(1)*q(1) - r(2)*q(2) - r(3)*q(3) - r(4)*q(4);
qr(2) = r(1)*q(2) + r(2)*q(1) - r(3)*q(4) + r(4)*q(3);
qr(3) = r(1)*q(3) + r(2)*q(4) + r(3)*q(1) - r(4)*q(2);
qr(4) = r(1)*q(4) - r(2)*q(3) + r(3)*q(2) + r(4)*q(1);
end