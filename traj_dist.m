function [xc, yc, xin, yin, xout, yout] = traj_dist(xc, yc, xin, yin, xout, yout, Ss, fig, tol)
    % This function returns a trajectory which has the same shape as the
    % input trajectory, but the samples have the same constant distance Ss.
    % The user can set also an allowed tolerance for the distance between
    % the samples, through the optional argument tol, whose default value is 1e-6.
    % 
    % INPUTS
    % (xc, yc)          samples of the center line
    % (xin, yin)        samples of the inner border
    % (xout, yout)      samples of the outer border
    % Ss                required distance between the samples
    % fig               a string to specify if the user want the plot of the trajectories (both input and output)
    % tol               (optional) allowed tolerance for the distance
    % 
    % OUTPUTS
    % (xc, yc)          samples of the center line
    % (xin, yin)        samples of the inner border
    % (xout, yout)      samples of the outer border

    if ~(strcmp(fig,'y') || strcmp(fig,'n'))
        error('The argument fig must be a string with characters y or n. If you want to plot the figures of the trajectory, please set fig to y, otherwise set it to n');
    end
    % tol is an optional argument
    if nargin < 9
        t = 1e-6;
    else
        t = tol;
    end
    if strcmp(fig, 'y')
        figure
        plot(xc, yc, '-o')
        hold on
        plot(xin, yin, '-og')
        hold on
        plot(xout, yout, '-og')
        grid on
        axis equal
        title('input trajectory')
    end

    for i = 1:length(xc)-1
        dc(i) = norm([xc(i+1)-xc(i), yc(i+1)-yc(i)]);
    end
    for i = 1:length(xin)-1
        din(i) = norm([xin(i+1)-xin(i), yin(i+1)-yin(i)]);
    end
    for i = 1:length(xout)-1
        dout(i) = norm([xout(i+1)-xout(i), yout(i+1)-yout(i)]);
    end
    s_ref = min([dc, din, dout]);
    flag = 0;
    if Ss < s_ref
        % the sampling distance is ok => use only the "easy" interpolation        
        Sin = Ss;
    else
        % the output trajectory must be computed in two step: first the
        % function traj_dist returns a trajectory of "finer" points, and
        % then the function discretizer returns the "real" output trajectory
        i = 1;
        while flag < 1
            Sin = Ss/i;
            if Sin < s_ref
                flag = 1;
            else
                i = i+1;
            end
        end
    end
    
    [xc, yc, ~] = finer_interpolation(xc,yc,Sin,t);
    [xin, yin, ~] = finer_interpolation(xin,yin,Sin,t);
    [xout, yout, ~] = finer_interpolation(xout,yout,Sin,t);
    if strcmp(fig, 'y')
        figure
        plot(xc, yc, '-o')
        hold on
        plot(xin, yin, '-og')
        hold on
        plot(xout, yout, '-og')
        grid on
        axis equal
        if flag > 0
            title('finer trajectory (not the output)')
        else
            title('output trajectory')
        end
    end
    if flag > 0
        xc = output_interpolation(xc, Ss, Sin);
        yc = output_interpolation(yc, Ss, Sin);
        xin = output_interpolation(xin, Ss, Sin);
        yin = output_interpolation(yin, Ss, Sin);
        xout = output_interpolation(xout, Ss, Sin);
        yout = output_interpolation(yout, Ss, Sin);
        if strcmp(fig, 'y')
            figure
            plot(xc, yc, '-o')
            hold on
            plot(xin, yin, '-og')
            hold on
            plot(xout, yout, '-og')
            grid on
            axis equal
            title('output trajectory')
        end
    end
end

function [x, y, d] = finer_interpolation(xc, yc, Ts, toll)
    % This function produce a trajectory with a given distance between the sample points
    % 
    % INPUTS
    % (xc, yc)      input trajectory
    % Ts            distance between points
    % toll          allowed tolerance
    % 
    % OUTPUTS
    % (x, y)        output trajectory
    % d             vector of the distances between the points of the output trajectory
    % 
    % A fundamental assumption is that the distance between each point in
    % (xc, yc) is higher than Ts, for all the points in the input trajectory
    
    [x, y, spre] = rough_interpolation(xc,yc,Ts);           % first attempt
    for i = 1:length(x)-1
        d(i) = norm([x(i+1)-x(i), y(i+1)-y(i)]);
    end
    % The following cicle assures that the distance between the intersamples points is Ts,
    % with a given allowable tolerance
    while or(max(d) > Ts*(1+toll), min(d) < Ts*(1-toll))
        j = 1;                      % index for s, variable defined in the following lines
        s(1) = spre(1); 
                                    % s is the query array related to the trajectory with sampling distance Ts
                                    % spre is the "old" query array, which must be corrected
        for i = 1:length(d)
            if or(d(i) < Ts*(1-toll), d(i) > Ts*(1+toll))
                % The distance is not ok => the query point must change
                if abs(s(j)-spre(i)) <= Ts/d(i)*(spre(i+1)-spre(i))
                    % The scaling factor between s and spre is Ts/d(i): if Ts<d(i), it means that 
                    % the "old" query point spre(i+1) must be closer to the previous one spre(i),
                    % otherwise if Ts>d(i) spre(i+1) must be further than spre(i). This information must be stored in s(j+1)
                    s(j+1) = s(j) + Ts/d(i)*(spre(i+1)-spre(i));
                else
                    % The distance between s(j) and spre(i) is higher than the "distance" that 
                    % must exists between each point of s(j) (to have the prescribed distance Ts).
                    % In this case the index between s(j) and spre(i) should different, i.e. we should add a new point
                    % s(j+1), and we could give s(j+1) the value spre(i): the distance between s(j+1) and s(j)
                    % would not be exactly Ts, but would not be too different from it,
                    % and shifting the indexes assures that the real output trajectory has enough point 
                    % to have an inter-samples distance equal to Ts
                    s(j+1) = spre(i);
                end
            else
                % The distance is ok => no need to compute a new query point
                s(j+1) = s(j) + spre(i+1) - spre(i);
            end
            j = j+1;
        end
        % It may happen that the last point of s is "too far" from the end of the original trajectory (xc,yc),
        % and in that way during the while cycle last point produced by the interpolation will be much further
        % than the last point of the original trajectory, and overall the output trajectory will be smaller.
        % The following check avoid this issue
        if length(xc)-1-s(end) > s(end)-s(end-1)
            s(end+1) = s(end) + s(end)-s(end-1);
        end
        spre = s;
        x = interp1(0:length(xc)-1,xc,s,'spline');
        y = interp1(0:length(yc)-1,yc,s,'spline');
        for i = 1:length(x)-1
            d(i) = norm([x(i+1)-x(i), y(i+1)-y(i)]);
        end
    end
end

function [x, y, s] = rough_interpolation(xc, yc, Ts)
    % This function makes a first interpolation of a trajectory, with a
    % desired distance between the samples (even though this last constraint is no respected
    % so precisely, and to have a more robust constraint use instead the function traj_dist)
    % 
    % INPUTS
    % (x, y)        is the interpolated trajectory
    % s             are the query points at which the trajectory has been interpolated
    % 
    % OUTPUTS
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
function t_out = output_interpolation(t_in, Ss_new, Ss_old)
    % This function is called when in the input trajectory exists a pair of points p1 = (x1,y1) and p2 = (x2,y2)
    % such that the distance between p1 and p2 is lower than Ss, so that the function finer_interpolation cannot work.
    % To avoid this issue there's a first run of the function finer_interpolation, that provides a trajectory
    % with distance Sin, which is at the same time an integer multiple of Ss and is smaller than the lowest distance
    % of the input trajectories.
    % This funtion starts from the trajectory provided by finer_interpolation and selects the sample point such that 
    % the distance between them is Ss.
    % This function is a simple for cicle along one direction of the given trajectory.
    % 
    % INPUTS 
    % t_in              input trajectory (only one dimension)
    % Ss_new            real distance required by the user (Ss in traj_dist)
    % Ss_old            finer distance used by finer_interpolation
    % 
    % OUTPUTS
    % t_out            output trajectory (only one dimension) 

    j = 1;
    for i = 1:Ss_new/Ss_old:length(t_in)
        t_out(j) = t_in(i);
        j = j+1;
    end
    
end
