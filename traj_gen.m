function [chi, k, s] = traj_gen(x, y, step)
    % chi       is the angle tangent to the trajectory
    % k         is the curvature (the inverse of the curvature radius)
    % s         is the real step size, that will be as close as possible to the given one
    % 
    % (x, y)    is the coordinate of the trajectory
    % step      is the given step size

    [~, ~, x_temp,y_temp,~] = ascissa_curvilinea(x,y,step);
                    % ascissa_curvilinea generates the waypoints on the trajectory at a given step
    x_temp = smooth(x_temp);
    y_temp = smooth(y_temp);
    x_true = x_temp';
    y_true = y_temp';
    clear x_temp y_temp
    chi = zeros(size(x_true));
    k = zeros(size(x_true));
    [chi(2:end), k(2:end),s] = traj_param(x_true,y_true);
end