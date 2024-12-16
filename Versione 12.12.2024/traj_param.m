function [chi, k, s] = traj_param(x, y)
    % This function provides the parameters of a trajectory starting from
    % the coordinates of the trajectory
    % 
    % chi       is the angle tangent to the trajectory
    % k         is the curvature (the inverse of the curvature radius)
    % s         is the real step size, that will be as close as possible to the given one
    % 
    % (x, y)    is the coordinate of the trajectory
    
    dy = diff(y);
    dx = diff(x);
    chi = atan2(dy,dx);                 % tangent angle of the trajectory wrt the global coordinates
    alpha(1) = 0;
    k(1) = 0;
    for i = 2:length(x)-1
        [b, ~] = cart2pol(x(i+1)-x(i),y(i+1)-y(i));
        [c, ~] = cart2pol(x(i+1)-x(i-1),y(i+1)-y(i-1));
        [~, a] = cart2pol(x(i)-x(i-1),y(i)-y(i-1));
        alpha(i) = b-c;
        k(i) = 2*sin(alpha(i))/a;       % curvature, computed using the geometrical properties of the trajectory
    end
    
    for i=2:length(x)
        p1 = [x(i),y(i)];
        p2 = [x(i-1),y(i-1)];
        s(i) = norm(p2-p1);             % the actual step size of the input trajectory
    end
end