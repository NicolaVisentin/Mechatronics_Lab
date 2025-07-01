clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   NOTATION                                                              %
%       - zx, zy:       [x xp]', [y yp]', (states)                        %
%       - u_x, u_y:     alphax,alphay (control)                           %
%       - P:            final state weighting matrix                      %
%       - PP:           Riccati matrix                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data and parameters definition

% Data of the mechanical system

m=0.264;            % mass of the ball
r=0.02;             % radius of the ball
g=9.81;             % gravity
J=2/5*m*r^2;        % inertia of the ball
b=0.0245;           % arm of the motor
lx=0.182;           % longer edge HALF length
ly=0.150;           % shorter edge HALF length

Cx=m*g*b/lx/(m+J/r^2);
Cy=m*g*b/ly/(m+J/r^2);

data_sys=[m r g J b lx ly];

% Load calibration

load('calibration.mat')

% LQR parameters (weights on delta_x and delta_u; final target)

p_delta_x = 0.1/(1)^2;
p_delta_y = 0.1/(1)^2;
p_delta_xp = 0;
p_delta_yp = 0;
r_delta_alphax = 0.05/(1)^2;
r_delta_alphay = 0.05/(1)^2;
q_delta_x = 0.7/(1)^2;
q_delta_y = 0.7/(1)^2;
q_delta_xp = 0;
q_delta_yp = 0;

Px_lqr=[p_delta_x 0; 0 p_delta_xp];     % Phi=0.5*(x(tf)-xf)'*P*(x(tf)-xf)
Rx_lqr=r_delta_alphax;                  % L=0.5*(u'*R*u+x'*Q*x)
Qx_lqr=[q_delta_x 0; 0 q_delta_xp];     % L=0.5*(u'*R*u+x'*Q*x)

Py_lqr=[p_delta_y 0; 0 p_delta_yp];     % Phi=0.5*(x(tf)-xf)'*P*(x(tf)-xf)
Ry_lqr=r_delta_alphax;                  % L=0.5*(u'*R*u+x'*Q*x)
Qy_lqr=[q_delta_y 0; 0 q_delta_yp];     % L=0.5*(u'*R*u+x'*Q*x)

%% Nominal trajectory

% load nominal trajectory

% load('u_z_zp_spiral.mat');   % contains: t_k, ux_nom, uy_nom, x_nom, y_nom, xp_nom, yp_nom

load('circle_traj.mat')
t_k=circle.t';
ux_nom=circle.ux;
uy_nom=circle.uy;
x_nom=circle.x;
y_nom=circle.y;
xp_nom=circle.xp;
yp_nom=circle.yp;

% load('track_traj.mat')
% t_k=track.t';
% ux_nom=track.ux;
% uy_nom=track.uy;
% x_nom=track.x;
% y_nom=track.y;
% xp_nom=track.xp;
% yp_nom=track.yp;

% plot nominal trajectory

figure

subplot(3,1,1)
hold on
plot(t_k,1000*x_nom,'b')
plot(t_k,1000*y_nom,'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('pos [mm]','Interpreter','latex')
title('Nominal x and y')
legend('$x^*$','$y^*$','Interpreter','LaTex','Location','best')
hold off

subplot(3,1,2)
hold on
plot(t_k,xp_nom,'b')
plot(t_k,yp_nom,'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('vel [m/s]','Interpreter','latex')
title('Nominal xp and yp')
legend('$\dot x^*$','$\dot y^*$','Interpreter','LaTex','Location','best')
hold off

subplot(3,1,3)
hold on
plot(t_k,rad2deg(ux_nom),'b')
plot(t_k,rad2deg(uy_nom),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
title('Nominal controls')
legend('$\alpha_x^*$','$\alpha_y^*$','Interpreter','LaTex','Location','best')
hold off

figure
hold on
plot(x_nom,y_nom)
rectangle('position',[-lx -ly 2*lx 2*ly],'LineWidth',2)
grid on
box on
xlabel('x')
ylabel('y')
title('Nominal trajectory')
axis([-1.2*lx 1.2*lx -1.2*ly 1.2*ly])
hold off

% re-sample nominal trajectory

% N=101;    % number of time intervals for re-sampling
% time=linspace(0,t_k(end),N);
time=t_k;
N=length(time);

xk=interp1(t_k,x_nom,time)';
yk=interp1(t_k,y_nom,time)';
xpk=interp1(t_k,xp_nom,time)';
ypk=interp1(t_k,yp_nom,time)';

zx=[xk xpk];
zy=[yk ypk];

ux=interp1(t_k,ux_nom,time)';
uy=interp1(t_k,uy_nom,time)';

Nz=2;
Nu=1;

%% LQR (for x direction)

% solve Riccati equation to find matrix PP.
%
%   ! backward integration
%   ! Riccati equation is a matrix equation, but ode45 only works with
%     vectors; we need to "unwrap" all matrices in input to ode45

options=odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,Nz^2));
t0=time(1);
tf=time(end);

PP_tf=Px_lqr;                % PP(tf) = d^2/dx^2(Phi)|tf = P
PP_tf_vect=PP_tf(1:end)';   % convert initial condition matrix into vector to input ode45

[t_PP,PP_vect]=ode45(@(t,PP) DRE(t,PP,Qx_lqr,Rx_lqr,zx,time,ux,time,Cx),[tf t0],PP_tf_vect,options);

% note that PP_vect is actually a (N x Nx^2) matrix, where the i-th row
% corresponds to a vector that represents the "unwrapped" Riccati matrix in
% the i-th time instant. Let's flip, re-sample and reshape PP_vect to
% obtain Riccati matrix time history in a (Nx x Nx x N) 3D matrix where each
% one of the N "slices" corresponds to the Nx x Nx Riccati matrix at a
% certain time instant

PP_vect=flipud(PP_vect);              % flip PP_vector (was backward integrated)
t_PP=flipud(t_PP);                    % flip t_PP vector (was backward integrated)
PP_vect=interp1(t_PP,PP_vect,time);   % re-sample PP_vect
PP=reshape(PP_vect.',Nz,Nz,[]);       % reshape PP_vect into PP as said before
PP=permute(PP,[2, 1, 3]);

% compute the gain matrix K time history ( K is a Nx by Nu by N 3D matrix; 
% each "slice" Nu x nx is the K matrix at one of the N time instants)

K=zeros(Nz,Nu,N);
B=zeros(Nz,Nu,N);
for ii=1:N
    B(:,:,ii)=dfdu(ux(ii),Cx);
    K(:,:,ii)=inv(Rx_lqr)*B(:,:,ii)'*PP(:,:,ii);     % K(t)=R^-1*B'(t)*P(t)
end

Kx_lqr=K;

%% LQR (for y direction)

% solve Riccati equation to find matrix PP.
%
%   ! backward integration
%   ! Riccati equation is a matrix equation, but ode45 only works with
%     vectors; we need to "unwrap" all matrices in input to ode45

options=odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,Nz^2));
t0=time(1);
tf=time(end);

PP_tf=Py_lqr;                % PP(tf) = d^2/dx^2(Phi)|tf = P
PP_tf_vect=PP_tf(1:end)';   % convert initial condition matrix into vector to input ode45

[t_PP,PP_vect]=ode45(@(t,PP) DRE(t,PP,Qy_lqr,Ry_lqr,zy,time,uy,time,Cy),[tf t0],PP_tf_vect,options);

% note that PP_vect is actually a (N x Nx^2) matrix, where the i-th row
% corresponds to a vector that represents the "unwrapped" Riccati matrix in
% the i-th time instant. Let's flip, re-sample and reshape PP_vect to
% obtain Riccati matrix time history in a (Nx x Nx x N) 3D matrix where each
% one of the N "slices" corresponds to the Nx x Nx Riccati matrix at a
% certain time instant

PP_vect=flipud(PP_vect);              % flip PP_vector (was backward integrated)
t_PP=flipud(t_PP);                    % flip t_PP vector (was backward integrated)
PP_vect=interp1(t_PP,PP_vect,time);   % re-sample PP_vect
PP=reshape(PP_vect.',Nz,Nz,[]);       % reshape PP_vect into PP as said before
PP=permute(PP,[2, 1, 3]);

% compute the gain matrix K time history ( K is a Nx by Nu by N 3D matrix; 
% each "slice" Nu x nx is the K matrix at one of the N time instants)

K=zeros(Nz,Nu,N);
B=zeros(Nz,Nu,N);
for ii=1:N
    B(:,:,ii)=dfdu(uy(ii),Cy);
    K(:,:,ii)=inv(Ry_lqr)*B(:,:,ii)'*PP(:,:,ii);     % K(t)=R^-1*B'(t)*P(t)
end

Ky_lqr=K;

%% Kalman filter design

% A, B and C matrices

A=[0 1; 0 0];
Bx=[0; Cx];
By=[0; Cy];
C=[1 0];

% Covariance matrices

Qxk=diag([100, 1000]);         % disturbance covariance
Qyk=diag([100, 1000]);         % disturbance covariance
Rxk=.1;          % noise covariance
Ryk=.1;          % noise covariance

% observer gain matrix

[~,~,Kx_obs]=care(A',C',Qxk,Rxk);
Kx_obs=Kx_obs';

[~,~,Ky_obs]=care(A',C',Qyk,Ryk);
Ky_obs=Ky_obs';

% observer state matrices

A_obs=A;
Bx_obs=[Bx Kx_obs];
By_obs=[By Ky_obs];
C_obs=eye(2);
D_obs=zeros(2);

%% Simulink

Kx_lqr=squeeze(Kx_lqr)';   % K(t) gain matrix
Ky_lqr=squeeze(Ky_lqr)';   % K(t) gain matrix
time=time';

ICx_obs=[0 0]';    % observer initial conditions (on delta_zx)
ICy_obs=[0 0]';    % observer initial conditions (on delta_zy)