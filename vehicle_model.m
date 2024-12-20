clear variables
close all
clc

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

% initialization
    % global time-based model
X = 0;
Y = 0;
Vx = 50/3.6;
beta = 0;               % side slip angle
psi = 0;                % yaw angle
r = 0;                  % yaw rate
% local space-based model
t = 0;
n = 0;
eta = 0;                % psi - chi (yaw angle - the tangent to the trajectory

z0 = [X;Y;Vx;beta;psi;r];
z0_rel = [t;n;Vx;beta;eta;r];
th = [m;Jz;a;b;Cf;Cr;rw;mu;Af;Al;Cx;Td_max;Td_min;C_R;rho;0;0;Cz];
    % the parameters take into account also the trajectory parameters, that for the time-based model are useless

vx_min = 0.5/3.6;    % the minimum velocity at whihc the simulation have to stop
%% input
% to be changed if we want to change the simulatede trajectory
    % we assume that the input are constant, so that it's easier to provide
    % the same input to both the models
Td = Td_max/50;
delta = 2*pi/180;

u = [Td;delta];
%% Simulation with RK2
Ts = 1e-3;
Tend = 500;
tvec = 0:Ts:Tend;
N = length(tvec);
zout = zeros(6,N);
uout = zeros(2,N);
zout(:,1) = z0;
uout(:,1) = [Td,delta];

tic 
for ind=2:N
    zprime = zout(:,ind-1)+Ts/2*vehicle_original(0,zout(:,ind-1),uout(:,ind-1),th);
    zout(:,ind) = zout(:,ind-1)+Ts*vehicle_original(0,zprime,uout(:,ind-1),th);
    uout(:,ind)= uout(:,1);
end
t_RK2=toc;

%% ideal trajectory generation
Ss = 1e-2;                  % sampling frequency in space domain
% the following "while" cycle should be activated if the input Td = 0.
% otherwise, the "while" cycle should be commented

% i = 1;
% while(zout(3,i) > vx_min)
%     i = i+1;
% end
i = N;

x = zout(1,1:i-1);
y = zout(2,1:i-1);
g = zeros(1,i);             % yaw angle from kinematic model
x_traj = zeros(1,i);
y_traj = zeros(1,i);
for j = 2:i
    g(j) = g(j-1)+zout(3,j-1)/(a+b)*tan(delta)*Ts;
    x_traj(j) = x_traj(j-1)+zout(3,j-1)*cos(g(j-1))*Ts;
    y_traj(j) = y_traj(j-1)+zout(3,j-1)*sin(g(j-1))*Ts;
end

[chi_true(2:end), k_true(2:end),s_true] = traj_gen(x_traj,y_traj,Ss);
% the above parameters are essential for the local model
%% ideal trajectory analysis
% figure
% plot(x,y)
% hold on
% plot(x_temp,y_temp)
% hold on
% plot(x_true,y_true)
% legend('time','spaced','smoothed')
% axis([-5e-3,20e-3,-1e-6,5e-6])
% title('constant sampled trajectory')
% 
% figure
% plot(chi)
% hold on
% plot(chi_temp)
% hold on
% plot(chi_true)
% legend('time','spaced','smoothed')
% title('chi')
% 
% figure
% plot(s)
% hold on
% plot(s_temp)
% hold on
% plot(s_true)
% legend('time','spaced','smoothed')
% title('s')
% 
% figure
% plot(k)
% hold on
% plot(k_temp)
% hold on
% plot(k_true)
% legend('time','spaced','smoothed')
% title('k')
%% local model simulation - RK2
zout_rel = zeros(6,length(s_true));
uout = zeros(2,length(s_true));
zout_rel(:,1) = [z0_rel];
uout(:,1) = [Td,delta];

for ind=2:(length(s_true)-2)
    % the parameters related to the trajectory change for time to time, so should be updated
    theta = [m;Jz;a;b;Cf;Cr;rw;mu;Af;Al;Cx;Td_max;Td_min;C_R;rho;chi_true(ind);k_true(ind);Cz];
    zprime_rel = zout_rel(:,ind-1)+s_true(ind)/2*vehicle(0,zout_rel(:,ind-1),uout(:,ind-1),theta);
    zout_rel(:,ind) = zout_rel(:,ind-1)+s_true(ind)*vehicle(0,zprime_rel,uout(:,ind-1),theta);
    uout(:,ind)= uout(:,1);
end 
%% computation of the XY trajectory obtained by the local model
xout_rel = zeros(size(chi_true));
yout_rel = zeros(size(chi_true));
% change of coordinates
for ind=2:length(chi_true)
    xout_rel(ind) = x_true(ind)-zout_rel(2,ind)*sin(chi_true(ind));
    yout_rel(ind) = y_true(ind)+zout_rel(2,ind)*cos(chi_true(ind));
end
%% plots
figure
plot(zout(3,:))
title('v_x in time')
figure
plot(zout_rel(3,:))
title('v_x in space')

% trajectories comparison
figure
plot(xout_rel,yout_rel),grid on, xlabel('X (m)'),ylabel('Y (m)')
hold on
plot(x,y)
hold on
plot(x_traj,y_traj)
hold on
plot(x_true,y_true)
legend('space','time','ideal', 'ideal smooth')
axis equal

figure
plot(zout_rel(1,:))
title('dt/ds')

figure
plot((yout_rel-y_true)./y_true)
title('relative variation of yout_{rel} wrt y_{true}')

figure
plot(xout_rel)
title('xout_{rel}')
figure
plot(yout_rel)
title('yout_{rel}')

figure
plot(zout_rel(2,:))
title('n in space')

figure
plot(zout_rel(5,:))
title('local yaw angle')