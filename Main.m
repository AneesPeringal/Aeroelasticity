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



epsilon = m/M;
% x_tilde = x / b;
% y = h/b; 
% x_alpha = S_alpha / (M * b);  % x_cg / b

 
M = [1         x_alpha       0 ;
    x_alpha   r_alpha^2      0;
    0            0           1 ];

C = [ zeta_h + epsilon*zeta+ beta* U_tilde , -epsilon*zeta*lambda                  , -epsilon*zeta         ;
    -nu*U_tilde - eta*zeta*lambda          ,  zeta_alpha+epsilon*zeta*lambda^2     , epsilon*zeta*lambda  ;
    -zeta                                  ,  zeta*lambda                         , zeta                 ];

K = [Omega^2 + epsilon*gamma      ,   beta*U_tilde^2-eta*gamma*lambda                 ,  -epsilon*gamma       ;
    -epsilon*gamma*lambda         , r_alpha^2-nu*U_tilde^2 + epsilon*gamma*lambda^2     ,  epsilon*gamma*lambda;
    -gamma                       , gamma*lambda                                       , gamma               ];

Fq = [xi_h*y^3+epsilon*xi(y-x_tilde-lambda*alpha)^3                ;
    xi_alpha*alpha^3+epsilon*lambda*xi(x_tilde+lambda*alpha-y)^3 ;
    xi*(x_tilde+lambda*alpha-y)^3                                ];
