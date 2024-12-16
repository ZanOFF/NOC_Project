function [zdot_rel]=vehicle(tau,z_rel,u,theta)
    % evolution  of the model
    % z_rel         are the space variables at the previous step
    % u             are the inputs
    % theta         are the parameters
    % 
    % zdot_rel      are the space variables computed at this istant
    %%
    % parameters
    g = 9.81;
    m = theta(1,1);
    Jz = theta(2,1);
    a = theta(3,1);
    b = theta(4,1);
    Cf = theta(5,1);             % stiffness coefficient - front
    Cr = theta(6,1);             % stiffness coefficient - rear
    rw = theta(7,1);             % tire radius
    mu = theta(8,1);
    Af = theta(9,1);             % front area of the vehicle
    Al = theta(10,1);            % lateral area of the vehicle
    Cx = theta(11,1);            % drag resistance
    C_R = theta(14,1);           % rolling resistance coefficient
    rho = theta(15,1);           % air density
    chi = theta(17,1);           % angle of the tangent to the trajectory with respect to the horizontal axis
    k = theta(18,1);             % dchi/ds = chi'
    Cz = theta(16,1);            % drag downforce 
    
    % states 
    t = z_rel(1,1);
    n = z_rel(2,1);             % normal to the trajectory
    Vx = z_rel(3,1);            % Vx(s)
    beta = z_rel(4,1);          % side slip angle as function of space, beta(s)
    eta = z_rel(5,1);           % psi - chi (yaw angle - the tangent to the trajectory)
    r = z_rel(6,1);             % yaw rate as function of space
    
    % inputs
    Td = u(1,1);
    delta = u(2,1);
    
    % disturbances
    W = 0;
    
    % all the quantities are computed using the time-based model, and to obtain
    % the space-based model we can simply multiply the timed-based model times
    % dt/ds, so for the chain rule of the derivatives we would obtain:
    % z(s)' = dz/dt * dt/ds = dz/ds
    
    % Kinematics
    Vy = tan(beta)*Vx;
    alphaf = atan((Vy+a*r)/Vx)-delta;
    alphar = atan((Vy-b*r)/Vx);
    zf = tan(alphaf);
    zr = tan(alphar);
    
    % Dynamics
    Fzf = m*g*b/(a+b)+1/2*Cz/2*Af*Vx^2*rho;
    Fzr = m*g*a/(a+b)+1/2*Cz/2*Af*Vx^2*rho;
    Fx = Td/rw;
    Fxd = 1/2*rho*Af*Cx*Vx^2;
    Fr_f = C_R*Fzf*Vx;
    Fr_r = C_R*Fzr*Vx;
    Fyd = 1/2*rho*Al*W^2;
    Fyf = min(mu*Fzf,max(-mu*Fzf, -Cf*zf+(Cf^2*abs(zf)*zf)/(3*mu*Fzf)-(Cf^3*zf^3)/(27*mu^2*Fzf^2)));
    Fyr = min(mu*Fzr,max(-mu*Fzr, -Cr*zr+(Cr^2*abs(zr)*zr)/(3*mu*Fzr)-(Cr^3*zr^3)/(27*mu^2*Fzr^2)));
    
    % sdot = (Vx*cos(eta)-Vy*sin(eta))/(1-n*k);
    
    zdot_rel(1,1) = (1-n*k)/(Vx*cos(eta)-Vy*sin(eta));                  % dt/ds
    tdot = zdot_rel(1,1);
            %  useful when passing from time-domain to space-domain
    zdot_rel(2,1) = (Vx*sin(eta)+Vy*cos(eta))*tdot;                     % n'
    zdot_rel(3,1) = (Fx-Fyf*sin(delta)-Fr_f-Fr_r-Fxd)/m*tdot;           % Vx'
    zdot_rel(4,1) = ((Fyf*cos(delta)+Fyr+Fyd)/(m*Vx)-r)*tdot;           % beta'
    zdot_rel(6,1) = (a*Fyf*cos(delta)-b*Fyr)/Jz*tdot;                   % r'
    zdot_rel(5,1) = r*tdot - k;                                         % eta'
end