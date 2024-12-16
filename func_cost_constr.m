function out = func_cost_constr(x, Ss, N, th1, th2, z0, n_max, n_min)

    u = x;
    
    s_int = Ss/5;
    
    z = zeros(6, N);
    z(:,1) = z0;

    for ind = 2:N
        % take the vehicle parameters th1 and the angle and curvature of
        % the ideal trajectory of the starting point
        theta = th2(:, ind-1);  % the assignment of the full vector each time to change 2 elements may be inefficient
        % compute the midpoint for RH2
        z(:,ind) = evolution(z(:,ind-1), Ss, [u(ind-1);u(ind-1+N)], s_int, th1, theta);
    end

    n_SIM = z(2,:);
    h = [-n_SIM'+n_max';
        n_SIM'+n_min'];

    % cost function
    % wt = 1;
    % wn = 1;
    % wvx = 0;
    % wb = 0;
    % we = 0;
    % wr = 0;
    % Pz = diag([wt, wn, wvx, wb, we, wr]);
    % wT = 0;
    % wd = 0;
    % Pu = diag([wT,wd]);
    % f = z'*Pz*z;
    f = z(1,end)^2+z(2,:)*z(2,:)';
    out = [f;h];
end