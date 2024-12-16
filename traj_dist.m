% function [x, y, d, c, it] = traj_dist(xc, yc, Ts, toll)
function [x, y, d] = traj_dist(xc, yc, Ts, toll)
    % This function produce a trajectory with a given distance between the
    % sample points
    % 
    % (xc, yc)      input trajectory
    % Ts            distance between points
    % toll          allowed tolerance
    % 
    % (x, y)        output trajectory
    % d             vector of the distances between the points of the output trajectory
    % 
    % A fundamental assumption is that the distance between each point in
    % (xc, yc) is higher than Ts, for all the point in the input trajectory
    
    [x, y, spre] = one_step_interpolation(xc,yc,Ts);
    % it = 0;
    for i = 1:length(x)-1
        d(i) = norm([x(i+1)-x(i), y(i+1)-y(i)]);
    end
    % the following cicle assures that the distance between the intersamples points is Ts,
    % with a given allowable tolerance
    while or(max(d) > Ts*(1+toll), min(d) < Ts*(1-toll))
        j = 1;          % index for s, variable defined in the following lines
        s(1) = spre(1); 
                        % s is the query array related to the trajectory with sampling distance Ts
                        % spre is the "old" query array, which must be corrected
        for i = 1:length(d)
            if or(d(i) < Ts*(1-toll), d(i) > Ts*(1+toll))
                % the distance is not ok => the query point must change
                if abs(s(j)-spre(i)) <= Ts/d(i)*(spre(i+1)-spre(i))
                    % the scaling factor between s and spre is Ts/d(i): if Ts<d(i), it means that 
                    % the "old" query point spre(i+1) must be closer to the previous one spre(i),
                    % otherwise if Ts>d(i) it's the opposite. This information must be stored in s(j+1)
                    s(j+1) = s(j) + Ts/d(i)*(spre(i+1)-spre(i));
                else
                    % the distance between s(j) and spre(i) is higher than the "distance" that 
                    % must exists between each point of s(j) (to have the prescribed distance Ts).
                    % This means that there's "too much space" between s(j) and spre(i), and 
                    % to cover this space 
                    s(j+1) = spre(i);
                end
            else
                % the distance is ok => no need to compute a new query point
                s(j+1) = s(j) + spre(i+1) - spre(i);
            end
            j = j+1;
        end
        if length(xc)-1-s(end) > s(end)-s(end-1)
            s(end+1) = s(end) + s(end)-s(end-1);
        end
        spre = s;
        % it = it+1;
        % c(it) = spre(end);
        x = interp1(0:length(xc)-1,xc,s,'spline');
        y = interp1(0:length(yc)-1,yc,s,'spline');
        for i = 1:length(x)-1
            d(i) = norm([x(i+1)-x(i), y(i+1)-y(i)]);
        end
    end
end

function [x, y, s] = one_step_interpolation(xc, yc, Ts)
    % This function makes a first interpolation of a trajectory, with a
    % desired distance between the samples (even though this last
    % constraint is no respected so precisely, and to have a more robust function use instead traj_dist)
    % 
    % (x, y)        is the interpolated trajectory
    % s             are the query points at which the trajectory has been interpolated
    % 
    % (xc, yc)      is the starting trajectory
    % Ts            is the intersample distance

    j = 1;
    % the following cicle computes the query points that will be used for the 'interp1' function 
    for i = 1:length(xc)-1 
        d = norm([xc(i+1)-xc(i), yc(i+1)-yc(i)]);
        temp = 0:1/floor(d/Ts):1;
        s_pre(j:j+length(temp)-1) = temp';
        j = length(s_pre)+1;
    end
    s_out = zeros(size(s_pre));
    
    % the following cicle makes the query points in a format useful for the 'interp1' function  
    i = 2;                          % index of the first attempt of the query point (the one from 0 to 1)
    for j = 2:length(s_out)         % index of the query points to give in input to interp1
        if i > length(s_pre)
            break
        end
        if s_pre(i) == 0
            i = i+1;
        end
        s_out(j) = s_out(j-1) + s_pre(i)-s_pre(i-1);
        i = i+1;
    end
    s = s_out(1:j-1);
    x = interp1(0:length(xc)-1,xc,s,'spline');
    y = interp1(0:length(yc)-1,yc,s,'spline');
end
