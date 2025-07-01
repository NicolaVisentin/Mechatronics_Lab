clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   NOTATION                                                              %
%       - zx,zy:        [x xp]', [y yp]' (states)                         %
%       - ux,uy:        alphax, alphay (control)                          %
%       - P:            final state weighting matrix                      %
%       - PP:           Riccati matrix                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data

lx=0.35/2;
ly=0.288/2;
b=0.025;
m=0.264;
g=9.81;
r=0.02;
J=2/5*m*r^2;

Cx=m*g*b/lx/(m+J/r^2);
Cy=m*g*b/ly/(m+J/r^2);

% load calibration

load('calibration.mat')

% infinite time LQR parameters (weights on delta_x and delta_u; final target)

r_delta_alphax = .01/(1)^2;
r_delta_alphay = .01/(1)^2;
q_delta_x = 1/(1)^2;
q_delta_y = 1/(1)^2;
q_delta_xp = 0;
q_delta_yp = 0;

Rx_lqr=r_delta_alphax;                  % L=0.5*(u'*R*u+x'*Q*x)
Qx_lqr=[q_delta_x 0; 0 q_delta_xp];     % L=0.5*(u'*R*u+x'*Q*x)

Ry_lqr=r_delta_alphay;                  % L=0.5*(u'*R*u+x'*Q*x)
Qy_lqr=[q_delta_y 0; 0 q_delta_yp];     % L=0.5*(u'*R*u+x'*Q*x)

% A, B and C matrices

A=[0 1; 0 0];
Bx=[0; Cx];
By=[0; Cy];
C=[1 0];

%%%%%%%%%% LQRI steady state %%%%%%%%%%

% Extended matrices (for steady-state error)

A_ex=[A, zeros(2,1); C, 0];
Bx_ex=[Bx; 0];
By_ex=[By; 0];

Qx_ex=[Qx_lqr, zeros(2,1); zeros(1,2), .0000001];
Qy_ex=[Qy_lqr, zeros(2,1); zeros(1,2), .0000001];

% design of the LQR using matlab function

[Kx_lqr,~,poles_conx]=lqr(A_ex,Bx_ex,Qx_ex,Rx_lqr);
[Ky_lqr,~,poles_cony]=lqr(A_ex,By_ex,Qy_ex,Ry_lqr);

%%%%%%%%%% Kalman filter design %%%%%%%%%%

% Covariance matrices

Qxk=diag([100, 1000]);    % disturbance covariance
Qyk=diag([100, 1000]);    % disturbance covariance
Rxk=.1;                   % noise covariance
Ryk=.1;                   % noise covariance

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

% observer initial conditions

ICx_obs=[0 0]'-[0.05 0]';    % observer initial conditions (on delta_zx=zx-zx_ref)
ICy_obs=[0 0]'-[0 0]';    % observer initial conditions (on delta_zy=zy-zy_ref)

%%%%%%%%%% Load trajectory for tracking %%%%%%%%%%

% Load trajectory for tracking

load('circle_traj.mat')
t_plan=circle.t;
x_plan=circle.x;
y_plan=circle.y;