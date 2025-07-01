function [PPp_vect]=DRE(t,PP_vect,Q,R,x,t_x,u,t_u,param)

% Models the Differential Riccati Equation to input in ode45:
%
%       [t_PP,PP_vect]=ode45(@(t,PP) DRE(t,PP,Q,R,xk,time,uk,time,param),[tf t0],PP_tf_vect,options);
%
% Outputs (of ode45):
%   - t_PP:        time vector where PP_vect is defined (! it's backward time)
%   - PP_vect:     P Riccati matrix (backward) time history in vector form: it's actually a 
%                  (N x Nx^2) matrix, where each i-th row corresponds to the P
%                  Riccati matrix "unwrapped" in a vector, at a given time instant 
% DRE.m function:
%
%       [PPp_vect]=DRE(t,PP_vect,Q,R,x,t_x,u,t_u,param)
%
% Inputs:
%   - Q:           weighting matrix on state time history
%   - R:           weighting matrix on control time history
%   - x:           nominal (optimal) states trajectory
%   - t_x:         time vector where x is defined
%   - u:           nominal (optimal) control time history
%   - t_u:         time vector where u is defined
%   - param:       parameters of the system (struct)
%   - t, PP_vect:  these inputs are handled by ode45 (see above)
% Outputs:
%   - PPp_vect:    this output is handled by ode45 (see above)
%
% ! remember that ode45 works with vectors, but DRE is a matrix equation: 
% we need to "unwrap" Riccati vector (PP_vect: input of ode45) into a matrix 
% (PP), solve the Riccati equation, and then converting the result (PPp) into
% a vector (PPp_vect: output of ode45).
% Notation:
%   - P:         P final state weighting matrix
%   - PP:        P(t) Riccati matrix
%   - PPp:       time derivative of P(t) Riccati matrix
%   - PP_vect:   P(t) Riccati matrix converted into vector
%   - PPp_vect:  first time derivative of P(t) Riccati matrix converted into
%                vector


% ode45 discretises ad cazzum: make sure to take x and u at the "correct" 
% time istant (coherent with PP which is being integrated)

xk=interp1(t_x,x,t);
uk=interp1(t_u,u,t);

% matrices A(tk) and B(tk)

A=dfdx();
B=dfdu(uk,param);

% convert PP_vect (ode45 input) into a matrix (to solve DRE)

Nx=length(Q);    % number of states
PP=zeros(Nx,Nx);
PP(1:end)=PP_vect(1:end);

% solve Riccati equation

PPp=-(A'*PP+PP*A-PP*B*inv(R)*B'*PP+Q);

% convert the resulting matrix into a vector

PPp_vect=PPp(1:end)';

end