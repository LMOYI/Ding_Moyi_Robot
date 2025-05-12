function R = rotation_matrix(n,theta)
w_hat = [0, -n(3), n(2);n(3), 0, -n(1); -n(2), n(1), 0];
R = eye(3) + w_hat*sind(theta) + w_hat*w_hat*(1-cosd(theta));
end