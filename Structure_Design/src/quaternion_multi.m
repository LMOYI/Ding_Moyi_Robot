function qr = quaternion_multi(q,r)
qr = [q(4)*r(1:3)+r(4)*q(1:3)+cross(q(1:3),r(1:3)); q(4)*r(4)-q(1:3)'*r(1:3)];

% qr = [q(1)*r(1)-q(2:4)*r(2:4)', q(1)*r(2:end)+r(1)*q(2:4)+cross(q(2:4),r(2:4))];

% if qr(4) == 0
%     qr = qr(1:3);
% end
% if qr(1:3)==zeros(3,1)
%     qr = qr(4);
% end
end