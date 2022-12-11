%  close all;
%  clear ,clc;
% 
% % inputs 
% x_alpha = 0.2;
% r_alpha = 0.5; 
% beta = 0.2; 
% nu = 0.08; 
% Omega = 0.5;
% zeta_alpha = 0.01; 
% zeta_h = 0.01;
% lambda =1;
% syms zeta gamma U_tilde
% 

% 
% % epsilon = m/M;
% epsilon = 0.05; % obtained from the paper
% % x_tilde = x / b;
% % y = h/b; 
% % x_alpha = S_alpha / (M * b);  % x_cg / b
% 
%  
% M = [1         x_alpha       0 ;
%     x_alpha   r_alpha^2      0;
%     0            0           1 ];
% 
% C = [ zeta_h + epsilon*zeta+ beta* U_tilde , -epsilon*zeta*lambda                  , -epsilon*zeta         ;
%     -nu*U_tilde - epsilon*zeta*lambda          ,  zeta_alpha+epsilon*zeta*lambda^2     , epsilon*zeta*lambda  ;
%     -zeta                                  ,  zeta*lambda                         , zeta                 ];
% 
% K = [Omega^2 + epsilon*gamma      ,   beta*U_tilde^2-epsilon*gamma*lambda                 ,  -epsilon*gamma       ;
%     -epsilon*gamma*lambda         , r_alpha^2-nu*U_tilde^2 + epsilon*gamma*lambda^2     ,  epsilon*gamma*lambda;
%     -gamma                       , gamma*lambda                                       , gamma               ];
% 
% Fq = [xi_h*y^3+epsilon*xi(y-x_tilde-lambda*alpha)^3                ;
%     xi_alpha*alpha^3+epsilon*lambda*xi(x_tilde+lambda*alpha-y)^3 ;
%     xi*(x_tilde+lambda*alpha-y)^3                                ];
% 
% %% Routh stability analysis
% %%
% syms z
% BB = M*z^2 +C*z +K;
% det(BB)
% 
% %%
% % We analyse the stability of the airfoil by itself. Therefore, we truncate
% % the matrices.
% 
% MM = M(1:2,1:2);
% KK = K(1:2,1:2);
% CC = C(1:2,1:2);
% 
% 
% B = [eye(2) eye(2);
%     -inv(MM)*KK -inv(MM)*CC];
% 
% syms S(s)
% S(s) = det(B-s*eye(size(B)))
% 
% 
% 
% %% symbolic analysis of the system
% syms x_alpha r_alpha beta nu Omega zeta_alpha zeta_h lambda epsilon zeta gamma U_tilde xi_h y alpha x_tilde xi xi_alpha
% 
% M = [1         x_alpha       0 ;
%     x_alpha   r_alpha^2      0;
%     0            0           1 ];
% 
% 
% C = [ zeta_h + epsilon*zeta+ beta* U_tilde , -epsilon*zeta*lambda                  , -epsilon*zeta         ;
%     -nu*U_tilde - epsilon*zeta*lambda          ,  zeta_alpha+epsilon*zeta*lambda^2     , epsilon*zeta*lambda  ;
%     -zeta                                  ,  zeta*lambda                         , zeta                 ];
% 
% K = [Omega^2 + epsilon*gamma      ,   beta*U_tilde^2-epsilon*gamma*lambda                 ,  -epsilon*gamma       ;
%     -epsilon*gamma*lambda         , r_alpha^2-nu*U_tilde^2 + epsilon*gamma*lambda^2     ,  epsilon*gamma*lambda;
%     -gamma                       , gamma*lambda                                       , gamma               ];
% 
% Fq = [xi_h*y^3+epsilon*xi*(y-x_tilde-lambda*alpha)^3                ;
%     xi_alpha*alpha^3+epsilon*lambda*xi*(x_tilde+lambda*alpha-y)^3 ;
%     xi*(x_tilde+lambda*alpha-y)^3                                ];
% 
% 
% syms z
% BB = M*z^2 +C^z +K;

%% Getting Routh

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
zeta = 0.15

syms gamma zeta U_tilde xi_h y x_tilde alpha xi xi_alpha z
% gamma_list = linspace(0.2,1,100);
% U_tilde_list = linspace(0.8,1.8,100);

epsilon = 0.05;
figure;
%  for zeta = 0.15:0.15
%     for j = 1:length(U_tilde_list)
%         gamma = gamma_list(i);
%         U_tilde = U_tilde_list(j);
        M = [1         x_alpha       0 ;
            x_alpha   r_alpha^2      0;
            0            0           1 ];
        
        C = [ zeta_h + epsilon*zeta+ beta* U_tilde , -epsilon*zeta*lambda                  , -epsilon*zeta         ;
            -nu*U_tilde - epsilon*zeta*lambda          ,  zeta_alpha+epsilon*zeta*lambda^2     , epsilon*zeta*lambda  ;
            -zeta                                  ,  zeta*lambda                         , zeta                 ];
        
        K = [Omega^2 + epsilon*gamma      ,   beta*U_tilde^2-epsilon*gamma*lambda                 ,  -epsilon*gamma       ;
            -epsilon*gamma*lambda         , r_alpha^2-nu*U_tilde^2 + epsilon*gamma*lambda^2     ,  epsilon*gamma*lambda;
            -gamma                       , gamma*lambda                                       , gamma               ];
        
        Fq = [xi_h*y^3+epsilon*xi*(y-x_tilde-lambda*alpha)^3                ;
            xi_alpha*alpha^3+epsilon*lambda*xi*(x_tilde+lambda*alpha-y)^3 ;
            xi*(x_tilde+lambda*alpha-y)^3                                ];
        
        BB = M*z^2 +C*z +K;
        charecteristic = collect(det(BB),z);
        Coefficients = coeffs(charecteristic,z);
        syms eps
        routh_matrix = routh(Coefficients,eps);
        routh_first_column = routh_matrix(:,1);
        fp = fimplicit(routh_first_column(6))
        xlim([0.8 1.8])
        ylim([0.2 1]);
        view([90 -90])
% end

%         if length(routh_first_column(routh_first_column>0))==7
%             A(i,j) = 1;
%             plot(gamma,U_tilde,'o','Color','black','MarkerFaceColor','k')
%             hold on
%         else
%             A(i,j) = 0;
%         end
%     end
% end



% %% plotting the region of stability
% for i = 1:6
%     fimplicit(routh_first_column(i))
%     hold on
% end
% xlim([0.8 1.8]);
% ylim([0.2 1]);


%% Plotting figure 3
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
%zeta = 0.15;
epsilon = 0.05;
syms zeta gamma U_tilde xi_h y x_tilde alpha xi xi_alpha z

M = [1         x_alpha       0 ;
    x_alpha   r_alpha^2      0;
    0            0           1 ];

C = [ zeta_h + epsilon*zeta+ beta* U_tilde , -epsilon*zeta*lambda                  , -epsilon*zeta         ;
    -nu*U_tilde - epsilon*zeta*lambda          ,  zeta_alpha+epsilon*zeta*lambda^2     , epsilon*zeta*lambda  ;
    -zeta                                  ,  zeta*lambda                         , zeta                 ];

K = [Omega^2 + epsilon*gamma      ,   beta*U_tilde^2-epsilon*gamma*lambda                 ,  -epsilon*gamma       ;
    -epsilon*gamma*lambda         , r_alpha^2-nu*U_tilde^2 + epsilon*gamma*lambda^2     ,  epsilon*gamma*lambda;
    -gamma                       , gamma*lambda                                       , gamma               ];

Fq = [xi_h*y^3+epsilon*xi*(y-x_tilde-lambda*alpha)^3                ;
    xi_alpha*alpha^3+epsilon*lambda*xi*(x_tilde+lambda*alpha-y)^3 ;
    xi*(x_tilde+lambda*alpha-y)^3                                ];

BB = M*z^2 +C*z +K;
charecteristic = collect(det(BB),z);
Coefficients = coeffs(charecteristic,z);
syms eps
routh_matrix = routh(Coefficients,eps);
routh_first_column = routh_matrix(:,1);
rfc = routh_first_column(end)
figure
gamma_list = 0.1:0.01:1.2;
zeta_list = 0.1:0.004:0.3;
for i = 1:length(gamma_list)
    for j = 1:length(zeta_list)
        func_temp = subs(rfc,[gamma zeta],[gamma_list(i), zeta_list(j)]);
        e5 = matlabFunction(func_temp,'Vars',U_tilde);
        U_tilde_found(i,j) = fzero(e5,1);
        imagesc('XData',gamma_list(i),'YData',zeta_list(j),'CData',U_tilde_found(i,j))
        hold on
    end
end
%%
% x_alpha = 0.2;
% r_alpha = 0.5;
% beta = 0.2;
% nu = 0.08;
% Omega = 0.5;
% zeta_alpha = 0.01;
% zeta_h = 0.01;
% epsilon =0.05;
% lambda = 1;
% syms zeta U_tilde gamma
% 
% 
% 
% a6 = r_alpha^2 - x_alpha^2;
% a5 = (zeta_h+(1+epsilon)*zeta)*r_alpha^2+(nu*x_alpha+beta*r_alpha^2)*U_tilde+zeta_alpha-x_alpha^2*zeta+(2*x_alpha+lambda)*epsilon*lambda*zeta;
% a4 = (Omega^2 +zeta*zeta_h+(1+epsilon)*gamma+1)*r_alpha^2 + (beta*zeta_alpha+(nu*x_alpha+r_alpha^2*beta)*zeta)*U_tilde - (beta*x_alpha+nu)*U_tilde^2 +(zeta_h +(1+epsilon)*zeta)*zeta_alpha+epsilon*lambda^2*zeta*zeta_h-gamma*x_alpha^2+(2*x_alpha+lambda)*epsilon*lambda*gamma;
% a3 = (zeta*Omega^2+(1+gamma)*zeta_h+(1+epsilon)*zeta)*r_alpha^2 -(nu*zeta_h+beta*x_alpha*zeta+(1+epsilon)*nu*zeta)*U_tilde^2+(beta*zeta*zeta_alpha+gamma*nu*x_alpha + r_alpha^2*beta*(1+gamma))*U_tilde+(zeta*zeta_alpha+epsilon*gamma*lambda^2)*zeta_h+(1+epsilon)*gamma*zeta_alpha + (zeta_alpha+epsilon*lambda^2*zeta)*Omega^2;
% a2 = ((1+gamma)*Omega^2+zeta*zeta_h+(1+epsilon)*gamma)*r_alpha^2-(nu*Omega^2+nu*zeta*zeta_h+ beta*gamma*x_alpha + (1+epsilon)*gamma*nu)*U_tilde^2 +beta*(gamma*zeta_alpha+r_alpha^2*zeta)*U_tilde+(zeta*zeta_alpha+epsilon*gamma*lambda^2)*Omega^2 + gamma*zeta_alpha*zeta_h;
% a1 = (zeta*Omega^2+gamma*zeta_h)*r_alpha^2-nu*(zeta*Omega^2+gamma*zeta_h)*U_tilde^2 + beta*gamma*r_alpha^2*U_tilde+ gamma*zeta_alpha*Omega^2;
% a0 = (r_alpha^2-nu*U_tilde^2)*gamma*Omega^2;
% 
% coeff = [a6 a5 a4 a3 a2 a1 a0];
% Routh = routh(coeff,eps);
% Routh_first_colomn = Routh(:,1);
% first_column_zeta = subs(Routh_first_colomn,zeta, 0.15);
% figure
% fimplicit(first_column_zeta(6));
%% figure 3
figure
clear U_tilde_found
gamma_list = 0:0.01:0.57;
zeta_list = 0.05:0.01:0.3;
for i = 1:length(gamma_list)
    parfor j = 1:length(zeta_list)
        func_temp = subs(routh_first_column,[gamma zeta],[gamma_list(i), zeta_list(j)]);
        e5 = matlabFunction(func_temp(6),'Vars',U_tilde);
        U_tilde_found(i,j) = fsolve(e5,0.97);
%         plot3(gamma_list(i),zeta_list(j),U_tilde_found(i,j),'.')
%         hold on
    end
end

[G,Z] = meshgrid(gamma_list,zeta_list); 
surf(G,Z,U_tilde_found')
colormap jet
shading interp
colorbar


%% Energy balance
for i = 1:length(gamma_list)
    for j = 1:length(zeta_list)
        plot3(gamma_list(i),zeta_list(j),U_tilde_found(i,j),'.')
        hold on
    end
end