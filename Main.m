 close all;
 clear ,clc;

% inputs 
x_alpha = 0.2;
r_alpha = 0.5; 
beta = 0.2; 
nu = 0.08; 
Omega = 0.5;
zeta_alpha = 0.01; 
zeta_h = 0.01;
lambda =1;
syms zeta gamma U_tilde



% epsilon = m/M;
epsilon = 0.05; % obtained from the paper
% x_tilde = x / b;
% y = h/b; 
% x_alpha = S_alpha / (M * b);  % x_cg / b

 
M = [1         x_alpha       0 ;
    x_alpha   r_alpha^2      0;
    0            0           1 ];

C = [ zeta_h + epsilon*zeta+ beta* U_tilde , -epsilon*zeta*lambda                  , -epsilon*zeta         ;
    -nu*U_tilde - epsilon*zeta*lambda          ,  zeta_alpha+epsilon*zeta*lambda^2     , epsilon*zeta*lambda  ;
    -zeta                                  ,  zeta*lambda                         , zeta                 ];

K = [Omega^2 + epsilon*gamma      ,   beta*U_tilde^2-epsilon*gamma*lambda                 ,  -epsilon*gamma       ;
    -epsilon*gamma*lambda         , r_alpha^2-nu*U_tilde^2 + epsilon*gamma*lambda^2     ,  epsilon*gamma*lambda;
    -gamma                       , gamma*lambda                                       , gamma               ];

Fq = [xi_h*y^3+epsilon*xi(y-x_tilde-lambda*alpha)^3                ;
    xi_alpha*alpha^3+epsilon*lambda*xi(x_tilde+lambda*alpha-y)^3 ;
    xi*(x_tilde+lambda*alpha-y)^3                                ];

%% Routh stability analysis

% We analyse the stability of the airfoil by itself. Therefore, we truncate
% the matrices.

MM = M(1:2,1:2);
KK = K(1:2,1:2);
CC = C(1:2,1:2);


B = [eye(2) eye(2);
    -inv(MM)*KK -inv(MM)*CC];

syms S(s)
S(s) = det(B-s*eye(size(B)))



