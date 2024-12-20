function [n_min, n_max, xtr, ytr] = bound_gen(xin, yin, xout, yout, xc, yc)
    % This function provides the physical limits of the distance from the
    % path to follow. These physical limits are computedd considering the
    % inner and outer boundaries for the trajectory, and the trajectory to
    % follow
    % 
    % (xin, yin)        inner border
    % (xout, yout)      outer border
    % (xc, xc)          center line
    % 
    % n_min             vector of the minimum allowable "n"
    % n_max             vector of the maximum allowable "n"

    %% selecting the start and the end of the useful trajectory
    [i_sin, t_sin] = bound_reg(xc, yc, xin, yin, 's','in');
    [i_ein, t_ein] = bound_reg(xc, yc, xin, yin, 'e','in');
    [i_sout, t_sout] = bound_reg(xc, yc, xout, yout, 's','out');
    [i_eout, t_eout] = bound_reg(xc, yc, xout, yout, 'e','out');
    
    if strcmp(t_sin, 'b')
        i_in = i_sin;
        if strcmp(t_sout, 'b')
            i_out = i_sout;
            i_c = 1;
        else
            i_out = 1;
            i_c = i_sout;
        end
    else
        i_in = 1;
        if strcmp(t_sout, 'b')
            i_out = i_sout;
            i_c = 1;
        else
            i_out = 1;
            i_c = max(i_sin, i_sout);
        end
    end

    if strcmp(t_ein, 'b')
        e_in = i_ein-1;
        if strcmp(t_eout, 'b')
            e_out = i_sout-1;
            e_c = 0;
        else
            e_out = 0;
            e_c = i_eout-1;
        end
    else
        e_in = 0;
        if strcmp(t_eout, 'b')
            e_out = i_eout-1;
            e_c = 0;
        else
            e_out = 0;
            e_c = max(i_ein, i_eout) -1;
        end
    end

    xin = xin(i_in:end-e_in);
    yin = yin(i_in:end-e_in);
    xout = xout(i_out:end-e_out);
    yout = yout(i_out:end-e_out);
    xc = xc(i_c:end-e_c);
    yc = yc(i_c:end-e_c);
    xtr = xc;
    ytr = yc;
    %%
    [chi, ~, ~] = traj_param(xc,yc);
    chi = [0, chi];
    j = 1;
    % inner border
    for i = 1:length(xc)
        a_in = chi(i) + pi/2;
        if a_in < 0
            a_in = a_in + 2*pi;
        end
        [vl, ~] = cart2pol(xin(j)-xc(i), yin(j)-yc(i));
        if vl < 0
            vl = vl+2*pi;
        end
        [vr, ~] = cart2pol(xin(j+1)-xc(i), yin(j+1)-yc(i));
        if vr < 0
            vr = vr+2*pi;
        end
        if and(a_in-vl < 0, a_in-vr < 0)
            while and(a_in-vl < 0, a_in-vr < 0)
                j = j+1;
                if j > length(xin)-1
                    break
                end
                [vl, ~] = cart2pol(xin(j)-xc(i), yin(j)-yc(i));
                if vl < 0
                    vl = vl+2*pi;
                end
                [vr, ~] = cart2pol(xin(j+1)-xc(i), yin(j+1)-yc(i));
                if vr < 0
                    vr = vr+2*pi;
                end
            end
        end
        if j < length(xin)    
            s = sin(a_in);
            c = cos(a_in);
            n_max(i) = (xin(j)*c+yin(j)*s+xin(j+1)*c+yin(j+1)*s-2*(xc(i)*c+yc(i)*s))/2;
        else
            j = 1;
        end
    end
    % outer border

    for i=1:length(xc)
        a_out = chi(i) - pi/2- 2*pi;
        [vl, ~] = cart2pol(xout(j)-xc(i), yout(j)-yc(i));
        
        vl = vl-2*pi;
        
        [vr, ~] = cart2pol(xout(j+1)-xc(i), yout(j+1)-yc(i));
        
        vr = vr-2*pi;
        k=j;
        % check if the boundaries from j to the end satisfy the condition
        while ~and(a_out-vl > 0, a_out-vr <0)
            j = j+1;
            if j > length(xout)-1
                break
            end
            [vl, ~] = cart2pol(xout(j)-xc(i), yout(j)-yc(i));
            vl = vl-2*pi;
            [vr, ~] = cart2pol(xout(j+1)-xc(i), yout(j+1)-yc(i));
            vr = vr-2*pi;
        end
        %check if the bound. from j to the start satisfy the condition
        if(j>length(xout)-1)
            while ~and(a_out-vl > 0, a_out-vr <0)
                k = k-1;
                if k < 1
                    break
                end
                [vl, ~] = cart2pol(xout(k)-xc(i), yout(k)-yc(i));
                vl = vl-2*pi;
                [vr, ~] = cart2pol(xout(k+1)-xc(i), yout(k+1)-yc(i));
                vr = vr-2*pi;
            end
        end
        if and(j > length(xout)-1, k<1)
            break
        end
        s = sin(a_out);
        c = cos(a_out);
        if j>length(xout)-1
            j=k;
        end
            
        n_min(i) = (xout(j)*c+yout(j)*s+xout(j+1)*c+yout(j+1)*s-2*(xc(i)*c+yc(i)*s))/2;
    end

    % j = 1;
    % for i = 1:length(xc)
    %     a_out = chi(i) - pi/2- 2*pi;
    %     [vl, ~] = cart2pol(xout(j)-xc(i), yout(j)-yc(i));
    % 
    %     vl = vl-2*pi;
    % 
    %     [vr, ~] = cart2pol(xout(j+1)-xc(i), yout(j+1)-yc(i));
    % 
    %     vr = vr-2*pi;
    % 
    %     if ~and(a_out-vl > 0, a_out-vr <0)
    %         while ~and(a_out-vl > 0, a_out-vr < 0)
    %             j = j+1;
    %             if j > length(xout)-1
    %                 break
    % 
    %             end
    %             [vl, ~] = cart2pol(xout(j)-xc(i), yout(j)-yc(i));
    %             % if vl > 0
    %                 vl = vl-2*pi;
    %             % end
    %             [vr, ~] = cart2pol(xout(j+1)-xc(i), yout(j+1)-yc(i));
    %             % if vr > 0
    %                 vr = vr-2*pi;
    %             % end
    %         end
    %     end
    %     if j > length(xout)-1
    %         break
    %     end
    %     i
    %     j
    %     s = sin(a_out);
    %     c = cos(a_out);
    %     n_min(i) = (xout(j)*c+yout(j)*s+xout(j+1)*c+yout(j+1)*s-2*(xc(i)*c+yc(i)*s))/2;
    % end
end
