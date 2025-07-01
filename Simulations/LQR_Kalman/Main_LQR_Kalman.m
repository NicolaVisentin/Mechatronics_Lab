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

% LQR parameters (weights on delta_x and delta_u; final target)

p_delta_x = 0.01/(1)^2;
p_delta_y = 0.01/(1)^2;
p_delta_xp = 0.01;
p_delta_yp = 0.01;
r_delta_alphax = 0.01/(1)^2;
r_delta_alphay = 0.01/(1)^2;
q_delta_x = 1/(1)^2;
q_delta_y = 1/(1)^2;
q_delta_xp = 0.1;
q_delta_yp = 0.1;

Px_lqr=[p_delta_x 0; 0 p_delta_xp];     % Phi=0.5*(x(tf)-xf)'*P*(x(tf)-xf)
Rx_lqr=r_delta_alphax;                  % L=0.5*(u'*R*u+x'*Q*x)
Qx_lqr=[q_delta_x 0; 0 q_delta_xp];     % L=0.5*(u'*R*u+x'*Q*x)

Py_lqr=[p_delta_y 0; 0 p_delta_yp];     % Phi=0.5*(x(tf)-xf)'*P*(x(tf)-xf)
Ry_lqr=r_delta_alphax;                  % L=0.5*(u'*R*u+x'*Q*x)
Qy_lqr=[q_delta_y 0; 0 q_delta_yp];     % L=0.5*(u'*R*u+x'*Q*x)

%% Nominal trajectory

% load nominal trajectory

load('circle_traj.mat');   % contains: t, ux, uy, x, y, xp, yp

t_k=circle.t;
x_nom=circle.x;
y_nom=circle.y;
xp_nom=circle.xp;
yp_nom=circle.yp;
ux_nom=circle.ux;
uy_nom=circle.uy;

% load('track_traj.mat');   % contains: t, ux, uy, x, y, xp, yp
% 
% t_k=track.t;
% x_nom=track.x;
% y_nom=track.y;
% xp_nom=track.xp;
% yp_nom=track.yp;
% ux_nom=track.ux;
% uy_nom=track.uy;

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

% re-sample nominal trajectory

N=length(t_k);    % number of time intervals for re-sampling
time=linspace(0,t_k(end),N);

xk=interp1(t_k,x_nom,time)';
yk=interp1(t_k,y_nom,time)';
xpk=interp1(t_k,xp_nom,time)';
ypk=interp1(t_k,yp_nom,time)';

zxk=[xk xpk];
zyk=[yk ypk];

uxk=interp1(t_k,ux_nom,time)';
uyk=interp1(t_k,uy_nom,time)';

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

[t_PP,PP_vect]=ode45(@(t,PP) DRE(t,PP,Qx_lqr,Rx_lqr,zxk,time,uxk,time,Cx),[tf t0],PP_tf_vect,options);

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
    B(:,:,ii)=dfdu(uxk(ii),Cx);
    K(:,:,ii)=pagemtimes(inv(Rx_lqr)*B(:,:,ii)',PP(:,:,ii));     % K(t)=R^-1*B'(t)*P(t)
end

% compute stability matrix A(t) time history both for the controlled and
% the uncontrolled systems, then compute the poles in both cases

A_unc=zeros(Nz,Nz,N);
for ii=1:N
    A_unc(:,:,ii)=dfdx();
end
poles_unc=pageeig(A_unc);
poles_unc=squeeze(poles_unc)';

A_con=zeros(Nz,Nz,N);
for ii=1:N
    A_con(:,:,ii)=A_unc(:,:,ii)-B(:,:,ii)'*K(:,:,ii);    % A_con(t)=A_unc(t)-B'(t)*K(t)
end
poles_con=pageeig(A_con);
poles_con=squeeze(poles_con)';

% Plots and animations

animation_speed=10;

% plot the poles of the controlled and uncontrolled system

Fig1=figure;

subplot(2,2,1)
hold on
grid on
box on
xlabel('Re','Interpreter','latex')
ylabel('Im','Interpreter','latex')
title('Poles of the uncontrolled system')
axis([min(real(poles_unc(:,1)))-1, max(real(poles_unc(:,1)))+1, -min(abs(imag(poles_unc(:,1))))-1, max(abs(imag(poles_unc(:,1))))+1])

subplot(2,2,2)
hold on
grid on
box on
xlabel('Re','Interpreter','latex')
ylabel('Im','Interpreter','latex')
title('Poles of the controlled system')
axis([min(real(poles_con(:,1)))-1, max(real(poles_con(:,1)))+1, -max(abs(imag(poles_con(:,1))))-1, max(abs(imag(poles_con(:,1))))+1])

subplot(2,2,3)
hold on
plot(nan,nan,'b')
plot(nan,nan,'g')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('states','Interpreter','latex')
title('Nominal (optimal) trajectory - time')
legend('$x$','$\dot{x}$','Interpreter','latex','location','best')
axis([t0 tf min(min(zxk)) max(max(zxk))])

subplot(2,2,4)
hold on
grid on
box on
xlabel('$x$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
title('Nominal (optimal) trajectory - phase plane')
axis([min(zxk(:,1)) max(zxk(:,1)) min(zxk(:,2)) max(zxk(:,2))])

for ii=1:N/animation_speed
    
    clf(Fig1)
    
    subplot(2,2,1)
    hold on
    plot(real(poles_unc(ii*animation_speed,:)),imag(poles_unc(ii*animation_speed,:)),'kx')
    plot(real(poles_unc(1:ii*animation_speed,:)),imag(poles_unc(1:ii*animation_speed,:)),'k--')
    grid on
    box on
    xlabel('Re','Interpreter','latex')
    ylabel('Im','Interpreter','latex')
    title('Poles of the uncontrolled system')
    axis([min(real(poles_unc(:,1)))-1, max(real(poles_unc(:,1)))+1, -max(abs(imag(poles_unc(:,1))))-1, max(abs(imag(poles_unc(:,1))))+1])

    subplot(2,2,2)
    hold on
    plot(real(poles_con(ii*animation_speed,:)),imag(poles_con(ii*animation_speed,:)),'kx')
    plot(real(poles_con(1:ii*animation_speed,:)),imag(poles_con(1:ii*animation_speed,:)),'k--')
    grid on
    box on
    xlabel('Re','Interpreter','latex')
    ylabel('Im','Interpreter','latex')
    title('Poles of the controlled system')
    axis([min(real(poles_con(:,1)))-1, max(real(poles_con(:,1)))+1, -max(abs(imag(poles_con(:,1))))-1, max(abs(imag(poles_con(:,1))))+1])
    
    subplot(2,2,3)
    hold on
    plot(time(1:ii*animation_speed),zxk(1:ii*animation_speed,1),'b')
    plot(time(1:ii*animation_speed),zxk(1:ii*animation_speed,2),'g')
    grid on
    box on
    xlabel('t [s]','Interpreter','latex')
    ylabel('states','Interpreter','latex')
    title('Nominal (optimal) trajectory')
    legend('$x$','$\dot{x}$','Interpreter','latex','location','best')
    axis([t0 tf min(min(zxk)) max(max(zxk))])

    subplot(2,2,4)
    hold on
    plot(zxk(1:ii*animation_speed,1),zxk(1:ii*animation_speed,2))
    grid on
    box on
    xlabel('$x$','Interpreter','latex')
    ylabel('$\dot{x}$','Interpreter','latex')
    title('Nominal (optimal) trajectory - phase plane')
    axis([min(zxk(:,1)) max(zxk(:,1)) min(zxk(:,2)) max(zxk(:,2))])

    drawnow

end

% save gain matrix history

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

[t_PP,PP_vect]=ode45(@(t,PP) DRE(t,PP,Qy_lqr,Ry_lqr,zyk,time,uyk,time,Cy),[tf t0],PP_tf_vect,options);

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
    B(:,:,ii)=dfdu(uyk(ii),Cy);
    K(:,:,ii)=pagemtimes(inv(Ry_lqr)*B(:,:,ii)',PP(:,:,ii));     % K(t)=R^-1*B'(t)*P(t)
end

% compute stability matrix A(t) time history both for the controlled and
% the uncontrolled systems, then compute the poles in both cases

A_unc=zeros(Nz,Nz,N);
for ii=1:N
    A_unc(:,:,ii)=dfdx();
end
poles_unc=pageeig(A_unc);
poles_unc=squeeze(poles_unc)';

A_con=zeros(Nz,Nz,N);
for ii=1:N
    A_con(:,:,ii)=A_unc(:,:,ii)-B(:,:,ii)'*K(:,:,ii);    % A_con(t)=A_unc(t)-B'(t)*K(t)
end
poles_con=pageeig(A_con);
poles_con=squeeze(poles_con)';

% Plots and animations

animation_speed=500;

% plot the poles of the controlled and uncontrolled system

Fig1=figure;

subplot(2,2,1)
hold on
grid on
box on
xlabel('Re','Interpreter','latex')
ylabel('Im','Interpreter','latex')
title('Poles of the uncontrolled system')
axis([min(real(poles_unc(:,1)))-1, max(real(poles_unc(:,1)))+1, -min(abs(imag(poles_unc(:,1))))-1, max(abs(imag(poles_unc(:,1))))+1])

subplot(2,2,2)
hold on
grid on
box on
xlabel('Re','Interpreter','latex')
ylabel('Im','Interpreter','latex')
title('Poles of the controlled system')
axis([min(real(poles_con(:,1)))-1, max(real(poles_con(:,1)))+1, -max(abs(imag(poles_con(:,1))))-1, max(abs(imag(poles_con(:,1))))+1])

subplot(2,2,3)
hold on
plot(nan,nan,'b')
plot(nan,nan,'g')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('states','Interpreter','latex')
title('Nominal (optimal) trajectory - time')
legend('$y$','$\dot{y}$','Interpreter','latex','location','best')
axis([t0 tf min(min(zyk)) max(max(zyk))])

subplot(2,2,4)
hold on
grid on
box on
xlabel('$y$','Interpreter','latex')
ylabel('$\dot{y}$','Interpreter','latex')
title('Nominal (optimal) trajectory - phase plane')
axis([min(zyk(:,1)) max(zyk(:,1)) min(zyk(:,2)) max(zyk(:,2))])

for ii=1:N/animation_speed
    
    clf(Fig1)
    
    subplot(2,2,1)
    hold on
    plot(real(poles_unc(ii*animation_speed,:)),imag(poles_unc(ii*animation_speed,:)),'kx')
    plot(real(poles_unc(1:ii*animation_speed,:)),imag(poles_unc(1:ii*animation_speed,:)),'k--')
    grid on
    box on
    xlabel('Re','Interpreter','latex')
    ylabel('Im','Interpreter','latex')
    title('Poles of the uncontrolled system')
    axis([min(real(poles_unc(:,1)))-1, max(real(poles_unc(:,1)))+1, -max(abs(imag(poles_unc(:,1))))-1, max(abs(imag(poles_unc(:,1))))+1])

    subplot(2,2,2)
    hold on
    plot(real(poles_con(ii*animation_speed,:)),imag(poles_con(ii*animation_speed,:)),'kx')
    plot(real(poles_con(1:ii*animation_speed,:)),imag(poles_con(1:ii*animation_speed,:)),'k--')
    grid on
    box on
    xlabel('Re','Interpreter','latex')
    ylabel('Im','Interpreter','latex')
    title('Poles of the controlled system')
    axis([min(real(poles_con(:,1)))-1, max(real(poles_con(:,1)))+1, -max(abs(imag(poles_con(:,1))))-1, max(abs(imag(poles_con(:,1))))+1])
    
    subplot(2,2,3)
    hold on
    plot(time(1:ii*animation_speed),zyk(1:ii*animation_speed,1),'b')
    plot(time(1:ii*animation_speed),zyk(1:ii*animation_speed,2),'g')
    grid on
    box on
    xlabel('t [s]','Interpreter','latex')
    ylabel('states','Interpreter','latex')
    title('Nominal (optimal) trajectory')
    legend('$y$','$\dot{y}$','Interpreter','latex','location','best')
    axis([t0 tf min(min(zyk)) max(max(zyk))])

    subplot(2,2,4)
    hold on
    plot(zyk(1:ii*animation_speed,1),zyk(1:ii*animation_speed,2))
    grid on
    box on
    xlabel('$y$','Interpreter','latex')
    ylabel('$\dot{y}$','Interpreter','latex')
    title('Nominal (optimal) trajectory - phase plane')
    axis([min(zyk(:,1)) max(zyk(:,1)) min(zyk(:,2)) max(zyk(:,2))])

    drawnow

end

% save gain matrix history

Ky_lqr=K;

%% Kalman filter design

% A, B and C matrices

A=[0 1; 0 0];
Bx=[0; Cx];
By=[0; Cy];
C=[1 0];

% Covariance matrices

Qxk=diag([10, 100]);         % disturbance covariance
Qyk=diag([10, 100]);         % disturbance covariance
Rxk=.1;          % noise covariance
Ryk=.1;          % noise covariance

% observer gain matrix

[~,Kx_obs,polesx_kalm]=icare(A',C',Qxk,Rxk);
Kx_obs=Kx_obs';

[~,Ky_obs,polesy_kalm]=icare(A',C',Qyk,Ryk);
Ky_obs=Ky_obs';

% observer state matrices

A_obs=A;
Bx_obs=[Bx Kx_obs];
By_obs=[By Ky_obs];
C_obs=eye(2);
D_obs=zeros(2);

%% Simulink

% parameters

zx0=[0.05,0]';          % initial conditions
zy0=[0,0]';          % initial conditions

ICx_obs=[0 0]';    % observer initial conditions (on delta_zx)
ICy_obs=[0 0]';    % observer initial conditions (on delta_zy)

tf=30;

Kx_lqr=squeeze(Kx_lqr)';   % K(t) gain matrix
Ky_lqr=squeeze(Ky_lqr)';   % K(t) gain matrix
time=time';

% simulate and extract results

out=sim('LQR_Kalman.slx');

t_simul_temp=squeeze(out.tout);
x_simul=squeeze(out.x_sim);
alphax_simul=squeeze(out.alphax_sim);
x_ref=squeeze(out.x_ref);
alphax_ref=squeeze(out.alphax_ref);
y_simul=squeeze(out.y_sim);
alphay_simul=squeeze(out.alphay_sim);
y_ref=squeeze(out.y_ref);
alphay_ref=squeeze(out.alphay_ref);
delta_zx_hat=out.delta_zx_hat;
delta_zy_hat=out.delta_zy_hat;
xp_simul=squeeze(out.xp_sim);
yp_simul=squeeze(out.yp_sim);
x_meas=squeeze(out.x_meas);
y_meas=squeeze(out.y_meas);

t_simul=t_simul_temp(1):0.01:t_simul_temp(end);
x_simul=interp1(t_simul_temp,x_simul,t_simul)';
alphax_simul=interp1(t_simul_temp,alphax_simul,t_simul)';
x_ref=interp1(t_simul_temp,x_ref,t_simul)';
alphax_ref=interp1(t_simul_temp,alphax_ref,t_simul)';
y_simul=interp1(t_simul_temp,y_simul,t_simul)';
alphay_simul=interp1(t_simul_temp,alphay_simul,t_simul)';
y_ref=interp1(t_simul_temp,y_ref,t_simul)';
alphay_ref=interp1(t_simul_temp,alphay_ref,t_simul)';
delta_zx_hat=interp1(t_simul_temp,delta_zx_hat,t_simul);
delta_zy_hat=interp1(t_simul_temp,delta_zy_hat,t_simul);
xp_simul=interp1(t_simul_temp,xp_simul,t_simul)';
yp_simul=interp1(t_simul_temp,yp_simul,t_simul)';
x_meas=interp1(t_simul_temp,x_meas,t_simul)';
y_meas=interp1(t_simul_temp,y_meas,t_simul)';

xp_ref=interp1(t_k,xp_nom,t_simul)';
yp_ref=interp1(t_k,yp_nom,t_simul)';

%% Plots

% x, alpha_x and y, alpha_y in time

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*x_simul,'b')
plot(t_simul,1000*x_ref,'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [mm]','Interpreter','latex')
title('x position')
legend('$x$','$x_{ref}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,rad2deg(alphax_simul),'r')
plot(t_simul,rad2deg(alphax_ref),'r--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
title('Control in x')
legend('$\alpha_x$','$\alpha_{x,ref}$','Interpreter','latex','location','best')
hold off

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*y_simul,'b')
plot(t_simul,1000*y_ref,'b--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('y position')
legend('$y$','$y_{ref}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,rad2deg(alphay_simul),'r')
plot(t_simul,rad2deg(alphay_ref),'r--')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\alpha$ [deg]','Interpreter','latex')
title('Control in y')
legend('$\alpha_y$','$\alpha_{y,ref}$','Interpreter','latex','location','best')
hold off

% trajectory

figure
hold on
plot(1000*x_simul,1000*y_simul,'b')
plot(1000*x_simul(1), 1000*y_simul(1),'bo')
plot(1000*x_ref,1000*y_ref,'r')
rectangle('Position',[-1000*lx -1000*ly 1000*2*lx 1000*2*ly],'LineWidth',1.5)
grid on
box on
xlabel('x [mm]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('Trajectory of the ball')
legend('actual traj','','reference traj','Location','best')
axis([-1000*1.1*lx 1000*1.1*lx -1000*1.1*ly 1000*1.1*ly])
hold off

% observer

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*x_meas,'g')
plot(t_simul,1000*x_simul,'b')
plot(t_simul,1000*x_ref,'b--')
plot(t_simul,1000*(delta_zx_hat(:,1)+x_ref),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('x [mm]','Interpreter','latex')
title('x position')
legend('$x_{meas}$','$x$','$x_{ref}$','$\hat{x}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,xp_simul,'b')
plot(t_simul,xp_ref,'b--')
plot(t_simul,(delta_zx_hat(:,2)+xp_ref),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')
title('x velocity')
legend('$\dot{x}$','$\dot{x}_{ref}$','$\dot{\hat{x}}$','Interpreter','latex','location','best')
hold off

figure

subplot(2,1,1)
hold on
plot(t_simul,1000*y_meas,'g')
plot(t_simul,1000*y_simul,'b')
plot(t_simul,1000*y_ref,'b--')
plot(t_simul,1000*(delta_zy_hat(:,1)+y_ref),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('y [mm]','Interpreter','latex')
title('y position')
legend('$y_{meas}$','$y$','$y_{ref}$','$\hat{y}$','Interpreter','latex','location','best')
hold off

subplot(2,1,2)
hold on
plot(t_simul,yp_simul,'b')
plot(t_simul,yp_ref,'b--')
plot(t_simul,(delta_zy_hat(:,2)+yp_ref),'r')
grid on
box on
xlabel('t [s]','Interpreter','latex')
ylabel('$\dot{y}$ [m/s]','Interpreter','latex')
title('y velocity')
legend('$\dot{y}$','$\dot{y}_{ref}$','$\dot{\hat{y}}$','Interpreter','latex','location','best')
hold off

% animation

filename='animation.gif';
FigTag=figure;
for ii=1:100:length(t_simul)

    clf(FigTag)
    
    x_ii=x_simul(ii);
    y_ii=y_simul(ii);

    hold on
    grid on
    plot(1000*x_simul(1:ii),1000*y_simul(1:ii),'k--')
    rectangle('Position',[-1000*lx -1000*ly 1000*2*lx 1000*2*ly],'LineWidth',1.5);
    hp=scatter(1000*x_ii,1000*y_ii,'fill');
    hp.MarkerEdgeColor='k';
    hp.MarkerFaceColor='k';
    hp.LineWidth=4;
    xlabel('$x\;[mm]$','Interpreter','LaTex')
    ylabel('$y\;[mm]$','Interpreter','LaTex')
    title('Animation')
    axis([-1000*1.1*lx 1000*1.1*lx -1000*1.1*ly 1000*1.1*ly])
    box on
    drawnow

    % save animation
    frame=getframe(FigTag);  % [ottieni il frame corrente della figura]
    img=frame2im(frame);     % [converti il frame in immagine RGB]
    [imind,cm]=rgb2ind(img,256);  % [converti l'immagine in formato indicizzato con 256 colori]
    if ii==1
        % [crea il file GIF; imposta quante volte la vuoi loopare e il ritardo tra i loop]
        imwrite(imind, cm, filename, 'gif', 'Loopcount', 1, 'DelayTime', 1);
    else
        % [aggiungi i frame successivi al file GIF e imposta il frame rate]
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/30);
    end

end