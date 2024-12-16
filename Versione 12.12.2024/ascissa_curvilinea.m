function [chi, k_out, x_new,y_new,s] = ascissa_curvilinea(x,y,Ss)
    % chi               is the angle tangent to the trajectory
    % k_out             is the curvature (the inverse of the curvature radius)
    % s                 is the real step size, that will be as close as possible to the given one
    % (x_new, y_new)    is the new trajectory discretization, with step size s
    % 
    % (x, y)            is the coordinate of the trajectory
    % Ss                is the given step size

    [chi, k_out, s] = traj_param(x, y);
    % now we must change the step size into the prescibed one Ss
    x_new = zeros(1,length(x)*(round(max(s)/Ss)+1));
    y_new = zeros(1,length(x)*(round(max(s)/Ss)+1));
    x_new(1) = x(1);
    y_new(1) = y(1);
    i = 2;
    k = i;
    while i < length(x)
        if(s(i) > Ss)
            n = s(i);
            grad = [x(i),y(i)]-[x(i-1),y(i-1)];
            j = k;
            fix = n/round(n/Ss);
            for l =1:round(n/Ss)
                x_new(j) = x_new(j-1)+fix*grad(1)/norm(grad);
                y_new(j) = y_new(j-1)+fix*grad(2)/norm(grad);
                j=j+1;
            end
            k = j;
            x_new(k) = x(i);
            y_new(k) = y(i);
            i = i+1;
        else
            p2 = [x(i-1),y(i-1)];
            grad = [x(i),y(i)]-p2;
            while and(norm(grad) < Ss, i < length(x))
                i = i+1;
                grad = [x(i),y(i)]-p2;
            end
            x_new(k) = x(i);
            y_new(k) = y(i);
            k = k+1;
        end
    end
    x_new = x_new(1:k-1);
    y_new = y_new(1:k-1);
end