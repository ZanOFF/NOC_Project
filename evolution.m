function z_out = evolution(z0, s_end, u_in, s_int, th1, th2) 
    % u_in = 2x1
    % th2 = 2x1
    s = 0:s_int:s_end;
    th = [th1; th2];
    z = zeros(6,length(s));
    z(:,1) = z0;
    
    for ind = 2:length(s)
        % the parameters related to the trajectory change for time to time, so should be updated
        zprime = z(:,ind-1) + s_int/2*vehicle(0,z(:,ind-1),u_in,th);
        % compute the value of the state on the next step
        z(:,ind) = z(:,ind-1) + s_int*vehicle(0,zprime,u_in,th);
    end 
    z_out = z(:,end);
end