function [index, type] = bound_reg(xc, yc, xb, yb, pos, bound)
    % This function provides the index of the first useful value for the
    % overall computation. If the center line start too early, it provides
    % the first informative index of the center line. If the boundary is
    % too large, it provides the first informative index of the boundary.
    % Can evaluate at the start or at the end of the trajectories
    % 
    % (xc, yc)      trajectory of the center line
    % (xb, yb)      trajectory of the boundary
    % pos           position (start or end) of the trajectory: 's' or 'e'
    % bord          type of border: 'in' or 'out'
    % 
    % index
    % type          path for which the index is referring: 'c' (center line), 'b' (boundary)
    %               or 'n'/'e' (some errors has occoured) 

    [chi, ~, ~] = traj_param(xc,yc);
    if and(strcmp(pos,'s'), strcmp(bound,'in'))
        [index, type] = bound_s_in(xc, yc, xb, yb, chi);
    elseif and(strcmp(pos,'e'), strcmp(bound,'in'))
        [index, type] = bound_e_in(xc, yc, xb, yb, chi);
    elseif and(strcmp(pos,'s'), strcmp(bound,'out'))
        [index, type] = bound_s_out(xc, yc, xb, yb, chi);
    elseif and(strcmp(pos,'e'), strcmp(bound,'out'))
        [index, type] = bound_e_out(xc, yc, xb, yb, chi);
    else
        error('An error in the generation of the boundaries');
    end
end
% The structure of the following functions is the same
%   - compute the angles of the vector normal to the waypoint of the center line "a" and
%     the vectors that connect the waypoint with two adjecent points of the boundary
%   - select if the "moving" path is the center line or the boundary
%      - if the moving path is the center line, "a" have to be always normal to the waypoint, 
%        and the waypoint is cycling along the center path. Since the waypoint is changing, 
%        all the vectors have to be computed again
%      - if the moving path is the boundary, the waypoint is fixed and so "a" may not be computed again, 
%        while the other two vectors must be computed again with the new samples of the boundary

function [index, type] = bound_s_in(xc, yc, xb, yb, chi) 
    a = chi(1) + pi/2;
    [v1t, ~] = cart2pol(xb(1)-xc(1), yb(1)-yc(1));
    [v2t, ~] = cart2pol(xb(2)-xc(1), yb(2)-yc(1));
    if and(a-v1t > 0, a-v2t > 0)
        type = 'c';
        i = 2;
        while and(a-v1t > 0, a-v2t > 0)
            a = chi(i) + pi/2;
            [v1t, ~] = cart2pol(xb(1)-xc(i), yb(1)-yc(i));
            if v1t < 0
                v1t = v1t + 2*pi;
            end
            [v2t, ~] = cart2pol(xb(2)-xc(i), yb(2)-yc(i));
            if v2t < 0
                v2t = v2t + 2*pi;
            end
            i = i+1;
        end
        index = i-1;
    elseif and(a-v1t < 0, a-v2t < 0)
        type = 'b';
        i = 2;
        while and(a-v1t < 0, a-v2t < 0)
            [v1t, ~] = cart2pol(xb(i)-xc(1), yb(i)-yc(1));
            if v1t < 0
                v1t = v1t + 2*pi;
            end
            [v2t, ~] = cart2pol(xb(i+1)-xc(1), yb(i+1)-yc(1));
            if v2t < 0
                v2t = v2t + 2*pi;
            end
            i = i+1;
        end
        index = i-1;
    else
        type = 'n';
        index = 1;
        warning('Some issue may occurs in the boundary generation at the start of the inner');
    end
end

function [index, type] = bound_e_in(xc, yc, xb, yb, chi)
    a = chi(end) + pi/2;
    [v1t, ~] = cart2pol(xb(end)-xc(end), yb(end)-yc(end));
    [v2t, ~] = cart2pol(xb(end-1)-xc(end), yb(end-1)-yc(end));
    if and(a-v1t > 0, a-v2t > 0)
        type = 'b';
        i = 2;
        while and(a-v1t > 0, a-v2t > 0)
            v1t = v2t;
            if v1t < 0
                v1t = v1t + 2*pi;
            end
            [v2t, ~] = cart2pol(xb(end-i)-xc(end), yb(end-i)-yc(end));
            if v2t < 0
                v2t = v2t + 2*pi;
            end
            i = i+1;
        end
        index = i-1;
    elseif and(a-v1t < 0, a-v2t < 0)
        type = 'c';
        i = 1;
        while and(a-v1t < 0, a-v2t < 0)
            a = chi(end-i) + pi/2;
            [v1t, ~] = cart2pol(xb(end)-xc(end-i), yb(end)-yc(end-i));
            if v1t < 0
                v1t = v1t + 2*pi;
            end
            [v2t, ~] = cart2pol(xb(end-1)-xc(end-i), yb(end-1)-yc(end-i));
            if v2t < 0
                v2t = v2t + 2*pi;
            end
            i = i+1;
        end
        index = i-1;
    else
        type = 'n';
        index = 1;
        warning('Some issue may occurs in the boundary generation at the end of the inner');
    end
end

function [index, type] = bound_s_out(xc, yc, xb, yb, chi) 
    a = chi(1) - pi/2;
    [v1t, ~] = cart2pol(xb(1)-xc(1), yb(1)-yc(1));
    [v2t, ~] = cart2pol(xb(2)-xc(1), yb(2)-yc(1));
    if and(a-v1t < 0, a-v2t < 0)
        type = 'c';
        i = 2;
        while and(a-v1t < 0, a-v2t < 0)
            a = chi(i) + pi/2;
            [v1t, ~] = cart2pol(xb(1)-xc(i), yb(1)-yc(i));
            [v2t, ~] = cart2pol(xb(2)-xc(i), yb(2)-yc(i));
            i = i+1;
        end
        index = i-1;
    elseif and(a-v1t > 0, a-v2t > 0)
        type = 'b';
        i = 2;
        while and(a-v1t > 0, a-v2t > 0)
            [v1t, ~] = cart2pol(xb(i)-xc(1), yb(i)-yc(1));
            [v2t, ~] = cart2pol(xb(i+1)-xc(1), yb(i+1)-yc(1));
            i = i+1;
        end
        index = i-1;
    else
        type = 'n';
        index = 1;
        warning('Some issue may occurs in the boundary generation at the start of the outer');
    end
end

function [index, type] = bound_e_out(xc, yc, xb, yb, chi)
    a = chi(end) - pi/2;
    [v1t, ~] = cart2pol(xb(end)-xc(end), yb(end)-yc(end));
    [v2t, ~] = cart2pol(xb(end-1)-xc(end), yb(end-1)-yc(end));
    if and(a-v1t < 0, a-v2t < 0)
        type = 'b';
        i = 2;
        while and(a-v1t < 0, a-v2t < 0)
            v1t = v2t;
            [v2t, ~] = cart2pol(xb(end-i)-xc(end), yb(end-i)-yc(end));
            i = i+1;
        end
        index = i-1;
    elseif and(a-v1t > 0, a-v2t > 0)
        type = 'c';
        i = 1;
        while and(a-v1t > 0, a-v2t > 0)
            a = chi(end-i) + pi/2;
            [v1t, ~] = cart2pol(xb(end)-xc(end-i), yb(end)-yc(end-i));
            [v2t, ~] = cart2pol(xb(end-1)-xc(end-i), yb(end-1)-yc(end-i));
            i = i+1;
        end
        index = i-1;
    else
        type = 'n';
        index = 1;
        warning('Some issue may occurs in the boundary generation at the end of the outer');
    end
end