clear variables
close all
%% Model parameters
m = 798;
Jz = 400;
a = 3;
b = 2.63;
Cf = 56500;             % stiffness coefficient - front
Cr = 54500;             % stiffness coefficient - rear
rw = 0.2285;            % tire radius
mu = 1.7;
Af = 1.65;              % front area of the vehicle
Al = 3.75;              % lateral area of the vehicle
Cx = 0.7;               % drag resistance
C_R = 0.002;            % rolling resistance coefficient
rho = 1.2;              % air density
delta_max = 0.297;      % maximum sterring angle [rad]
Td_max = 950;           % maximum torque
Td_min = -3500;         % minimum torque
Cz = 2.5;               % drag downforce
th1 = [m;Jz;a;b;Cf;Cr;rw;mu;Af;Al;Cx;Td_max;Td_min;C_R;rho;Cz];
                        % model parameter in a vector (useful later)
%% Trajectory to follow and the boundaries 
% importing the coordinates
load('CL_traj.mat');
[xcl, ycl] = rimozione_outliers(coordinates);
clear coordinates
load('In_border.mat');
[xin, yin] = rimozione_outliers(coordinates);
clear coordinates
load('Out_border.mat');
[xout, yout] = rimozione_outliers(coordinates);
clear coordinates

% now the sample point of the trajectory have to have a prescribed distance
Ss = 0.5;                      % [m] - distance required btw samples
% for some reason works only for Ss=0.5
toll = 1e-3;                    % allowed tollerance for the distance
%[xc, yc, dc] = traj_dist(xcl,ycl,Ss,toll);
%[xin, yin, din] = traj_dist(xin,yin,Ss,toll);
%[xout, yout, dout] = traj_dist(xout,yout,Ss,toll);
[xc,yc,xin,yin,xout,yout]=traj_dist_new(xcl,ycl,xin,yin,xout,yout,Ss, 'y', toll);
%[xc,yc,xin,yin,xout,yout] = discretizer(xc,yc,xin,yin,xout,yout,2, Ss);
Send = length(xc);
% figure
% plot(dc),hold on,plot(dbin),hold on, plot(dbout)
% legend('cl','in','out')
% figure
% plot(xbin,ybin, xbout,ybout, xc,yc),grid on
% legend('inner','outer','ref')
% xlabel('X (m)'),ylabel('Y (m)')
% title('sequal distance')
% axis equal
% figure
% plot(xin,yin, xout,yout, xcl,ycl),grid on
% legend('inner','outer','ref')
% xlabel('X (m)'),ylabel('Y (m)')
% title('original')
% axis equal

[chi, k, s] = traj_param(xc, yc);

theta = [chi;k];                 % model parameter of the ideal trajectory
                                 % chi is angle of the tangent of the trajectory
                                 % k is the curvature for each discretization point
[n_min_tot, n_max_tot, xcl, ycl] = bound_gen(xin, yin, xout, yout, xc, yc);

%% Finite Horizon Optimal Control Problem
N = 80;                    % prediction horizon
s0 = 3;                     % fixed instants
% linear constraints
% x is u during all the prediction horizon
A = [eye(s0), zeros(s0, 2*N-s0);
    zeros(N-s0, 2*N);
    zeros(s0,N), eye(s0), zeros(s0,N-s0);
    zeros(N-s0, 2*N)];
b = zeros(2*N,1);
C = [-eye(2*N);
      eye(2*N)];

d = [-Td_max*ones(N,1);
     -delta_max*ones(N,1);
      Td_min*ones(N,1);
     -delta_max*ones(N,1)];

q = 2*N;
t0 = 0;
n0 = 1;
vx0 = 85/3.6;
beta0 = 0;
eta0 = 0;
r0 = 0;
Td0 = 0*eye(N,1);
delta0 = 0*eye(N,1);
myoptions               =   myoptimset;
myoptions.Hessmethod  	=	'BFGS';
myoptions.gradmethod  	=	'CD';
myoptions.graddx        =	2^-17;
myoptions.tolgrad    	=	1e-8;
myoptions.ls_beta       =	0.2;
myoptions.ls_c          =	0.1;
myoptions.ls_nitermax   =	1e2;
myoptions.nitermax      =	10;
myoptions.tolfun        =   1e-5;
myoptions.xsequence     =	'on';

z0 = [t0;n0;vx0;beta0;eta0;r0];
u0 = [Td0;delta0];
z = zeros(6,length(xcl));
z(:,1) = z0;
u_opt = zeros(2,length(xcl));
s_int = 5e-2;

for i = 1:length(xcl)-N
    xtr = xcl(i:N+i-1);
    ytr = ycl(i:N+i-1);
    % xbin = xin(i:N+i-1);
    % ybin = yin(i:N+i-1);
    % xbout = xout(i:N+i-1);
    % ybout = yout(i:N+i-1);
    th2 = theta(:, i:N+i-1);
    myoptions.outputfcn     =   @(x)Vehicle_traj(x,Ss,N,th1, th2,z(:,i), xtr, ytr,xin, yin, xout, yout);
                            % function handle, @(x) can be defined when calling the handel, while the other
                            % parameters must be defined when define the handle (i.e. here)
    n_min = n_min_tot(i:N+i-1);
    n_max = n_max_tot(i:N+i-1);
    tic
    [xstar,fxstar,iter,exitflag,xsequence] = myfmincon(@(x)func_cost_constr(x, Ss, N, th1, th2, z(:,i), n_max, n_min),u0,A,b,C,d,0,q,myoptions);
    t_exe(i) = toc
    u0 = [xstar(2:N); xstar(N); xstar(N+2:end); xstar(end)];            % we decided to fix the first input to have a softer u(t)
    
    b = [xstar(2:s0+1); zeros(N-s0,1);
        xstar(N+2:N+s0+1); zeros(N-s0,1)];
    u_opt(:,i) = [xstar(1); xstar(N+1)];
    z(:,i+1) = evolution(z(:,i), Ss, u_opt(:,i), s_int, th1, th2(:,1));
    close figure 1
end

for ind = 1:N-1
        xstar_tr(ind) = xcl(ind)-z(2,ind)*sin(chi(ind));
        ystar_tr(ind) = ycl(ind)+z(2,ind)*cos(chi(ind));
end

%%

figure
subplot(2,2,1),plot(xstar_tr,ystar_tr, xin,yin, xout,yout, xcl,ycl),grid on
legend('traj','inner','outer','ref')
xlabel('X (m)'),ylabel('Y (m)')
axis equal
subplot(2,2,2),plot(0:length(xcl)-1,z(1,:)),grid on
xlabel('Space (m)'),ylabel('Time to get this position (s)')
subplot(2,2,3),plot(0:length(xcl)-1,u_opt(2,:))
xlabel('Space (m)'),ylabel('Front steering angle (rad)')
subplot(2,2,4),plot(0:length(xcl)-1,u_opt(1,:))
xlabel('Space (m)'),ylabel('Driving torque (Nm)')

%%
figure
plot(0:length(xcl)-1,z(3,:)),grid on
xlabel('Space (s)'),ylabel('v_x (m/s)')
title('Longitudinal veloctiy in space')
