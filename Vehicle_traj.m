function [z] = Vehicle_traj(x,Ss,N,th1,th2, z0, xtr, ytr, xin, yin, xout, yout)
% Function that computes the trajectory of the vehicle model
% introduced in Lab. session A and returns the system state together with
% plots of the relevant quantities
    u = x;
        
    s_int = Ss/10;
    
    z = zeros(6,N);
    z(:,1) = z0;
    
    for ind = 2:N
        % the parameters related to the trajectory change for time to time, so should be updated
        theta = th2(:,ind-1);
        z(:,ind) = evolution(z(:,ind-1), Ss, [u(ind-1);u(ind-1+N)], s_int, th1, theta);
    end 
    
    chi = th2(1,:);
    for ind = 1:N-1
        xstar(ind) = xtr(ind)-z(2,ind)*sin(chi(ind));
        ystar(ind) = ytr(ind)+z(2,ind)*cos(chi(ind));
    end
    
    %% Plot (X,Y) trajectory and constraints
    s = 0:N-1;
    figure(1)
    subplot(2,3,1), plot(xstar,ystar, xin,yin, xout,yout, xtr,ytr),grid on
    %legend('traj','inner','outer','ref')
    xlabel('X (m)'),ylabel('Y (m)')
    axis equal
    subplot(2,3,2), plot(s,z(3,:)),grid on
    xlabel('Space (m)'),ylabel('Longitudinal speed (m/s)')
    subplot(2,3,3), plot(s,u(N+1:end)),grid on
    xlabel('Space (m)'),ylabel('Front steering angle (rad)')
    subplot(2,3,4), plot(s,u(1:N)),grid on
    xlabel('Space (m)'),ylabel('Driving torque (Nm)')
    subplot(2,3,5), plot(s,z(1,:)),grid on
    xlabel('Space (m)'),ylabel('Time to get this position (s)')
    subplot(2,3,6), plot(s,z(2,:)),grid on
    xlabel('Space (m)'),ylabel('n (m)')
end