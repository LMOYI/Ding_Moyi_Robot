
function  robot_plot(a, b, s, p, R, k, l)
%ROBOT_PLOT 给定结构参数a、b，驱动器长度s，位姿p、R，绘图模式k，移动方向l，绘制并联机构
% 输入参数
%       k：绘图模式参数，k=0，绘制基座，否则不绘制基座
if nargin <6
    l = [0,1,0]';
    k = 0;
elseif nargin == 6
    l = [0,1,0]';
end
p = p(1:3);
s = reshape(s,[1,6]);
% 创建图形窗口并设置视角
h=figure;
hold on;
grid on;
axis equal;
view(3);

% 计算移动后的基座铰点坐标
a_moved = a + s .* l; % 基座铰点沿l方向移动s距离

% % 绘制移动后基座（粉色半透明填充）
% fill3(a_moved(1,:), a_moved(2,:), a_moved(3,:), 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% 转换动平台铰点到世界坐标系
b_world = repmat(p, 1, 6) + R * b;

% 绘制动平台（黄色半透明填充）
fill3(b_world(1,:), b_world(2,:), b_world(3,:), 'y', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% 绘制移动后基座铰点（粉色小球）
scatter3(a_moved(1,:), a_moved(2,:), a_moved(3,:), 50, 'k', 'filled');

% 绘制动平台铰点（粉色小球）
scatter3(b_world(1,:), b_world(2,:), b_world(3,:), 50, 'k', 'filled');

% 绘制连杆（黑色实线）
for i = 1:6
    plot3([b_world(1,i), a_moved(1,i)], [b_world(2,i), a_moved(2,i)], [b_world(3,i), a_moved(3,i)], ...
        'k-', 'LineWidth', 2);
end

if k == 0
    % 绘制原基座（红色半透明填充）
    % fill3(a(1,:), a(2,:), a(3,:), 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    % 绘制原基座铰点（红色小球）
    scatter3(a(1,:), a(2,:), a(3,:), 50, 'r', 'filled');

    % 绘制移动前后基座铰点连线（红色实线）
    for i = 1:6
        plot3([a(1,i), a_moved(1,i)], [a(2,i), a_moved(2,i)], [a(3,i), a_moved(3,i)], ...
            'r-', 'LineWidth', 2);
    end


    % 绘制更新后的基座滑轨（灰色虚线）
    for i = 1:6
        start_point = a_moved(:,i) - l*0.25;
        end_point = start_point + l*0.5;
        plot3([start_point(1), end_point(1)], [start_point(2), end_point(2)], ...
            [start_point(3), end_point(3)], 'k--', 'LineWidth', 1);
    end
end
hold off;

% 设置坐标轴标签和标题
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Parallel Robot Mechanism with Motion');
end