close all
clear variables 

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

% parameters
Ts = 5e-1;              % distance required btw samples
toll = 1e-3;            % allowed tollerance for the distance
figure
plot(xcl, ycl, '-o')
hold on
plot(xin, yin, '-og')
hold on
plot(xout, yout, '-og')
grid on
axis equal
title('equal')

% [xc, yc, dc, cc, itc] = traj_dist(xcl,ycl,Ts,toll);
% [xin, yin, din, cin, itin] = traj_dist(xin,yin,Ts,toll);
% [xout, yout, dout, cout, itout] = traj_dist(xout,yout,Ts,toll);
[xc, yc, dc] = traj_dist(xcl,ycl,Ts,toll);
[xin, yin, din] = traj_dist(xin,yin,Ts,toll);
[xout, yout, dout] = traj_dist(xout,yout,Ts,toll);

figure
plot(xc, yc, '-o')
hold on
plot(xin, yin, '-og')
hold on
plot(xout, yout, '-og')
grid on
axis equal
title('distanced')
%%
[n_min, n_max, xcl, ycl] = bound_gen(xin, yin, xout, yout, xc, yc);
[chi, k, s] = traj_param(xcl, ycl);
chi = [0,chi];
x_upper = xcl + n_max.*cos(chi+pi/2*ones(size(chi)));
y_upper = ycl + n_max.*sin(chi+pi/2*ones(size(chi)));
x_lower = xcl + n_min.*cos(chi-pi/2*ones(size(chi)));
y_lower = ycl + n_min.*sin(chi-pi/2*ones(size(chi)));

figure
plot(xcl, ycl, '-o')
hold on
plot(xin, yin, '-og')
hold on
plot(xout, yout, '-og')
hold on
plot(x_upper, y_upper, ':xr')
hold on
plot(x_lower, y_lower, ':xr')
legend('cl','in','out')
grid on
axis equal

%%
[xc, yc, xin, yin, xout, yout] = discretizer(xc, yc, xin, yin, xout, yout, 2, Ts);
figure
plot(xc, yc, '-o')
hold on
plot(xin, yin, '-og')
hold on
plot(xout, yout, '-og')
grid on
axis equal
title('more space')
[n_min, n_max, xcl, ycl] = bound_gen(xin, yin, xout, yout, xc, yc);
[chi, k, s] = traj_param(xcl, ycl);
chi = [0,chi];
x_upper = xcl + n_max.*cos(chi+pi/2*ones(size(chi)));
y_upper = ycl + n_max.*sin(chi+pi/2*ones(size(chi)));
x_lower = xcl + n_min.*cos(chi-pi/2*ones(size(chi)));
y_lower = ycl + n_min.*sin(chi-pi/2*ones(size(chi)));
%%
figure
plot(xcl, ycl, '-o')
hold on
plot(xin, yin, '-og')
hold on
plot(xout, yout, '-og')
hold on
plot(x_upper, y_upper, ':xr')
hold on
plot(x_lower, y_lower, ':xr')
legend('cl','in','out')
grid on
axis equal
title('distanced')