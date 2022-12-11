%  close all;
 clear ,clc;
list = 0 : 0.25 : 1;
list = 0.5:0.5;
list = 1 * list ;
% list(end+1) = 0.462 ;
% for xx = list 


% inputs 
x_alpha = 0.2;
r_alpha = 0.5; 
beta = 0.2; 
nu = 0.08; 
Omega = 0.5;

zeta_alpha = 0.01; 
zeta_h = 0.01;
lambda =1;
% syms zeta gamma U_tilde
zeta = 0.11 ;
% zeta = xx ;
gamma = 0.462 ;
% gamma = xx ;
U_tilde =0.5* 1.255 ;
% U_tilde =1.3 ;

% epsilon = m/M;
epsilon = 0.05; % obtained from the paper
% epsilon = xx; % obtained from the paper
% x_tilde = x / b;
% y = h/b; 
% x_alpha = S_alpha / (M * b);  % x_cg / b

 xi_h = 0.06;  
 xi_alpha = 0.15;
 xi = 0.15;          % c/m*omega_alpha
%  xi = xx;          % c/m*omega_alpha


M = [1         x_alpha       0 ;
    x_alpha   r_alpha^2      0;
    0            0           1 ];

C = [ zeta_h + epsilon*zeta+ beta* U_tilde , -epsilon*zeta*lambda                  , -epsilon*zeta         ;
    -nu*U_tilde - epsilon*zeta*lambda          ,  zeta_alpha+epsilon*zeta*lambda^2     , epsilon*zeta*lambda  ;
    -zeta                                  ,  zeta*lambda                         , zeta                 ];

K = [Omega^2 + epsilon*gamma      ,   beta*U_tilde^2-epsilon*gamma*lambda                 ,  -epsilon*gamma       ;
    -epsilon*gamma*lambda         , r_alpha^2-nu*U_tilde^2 + epsilon*gamma*lambda^2     ,  epsilon*gamma*lambda;
    -gamma                       , gamma*lambda                                       , gamma               ];

% Fq = [xi_h*y^3+epsilon*xi(y-x_tilde-lambda*alpha)^3                ;
%     xi_alpha*alpha^3+epsilon*lambda*xi(x_tilde+lambda*alpha-y)^3 ;
%     xi*(x_tilde+lambda*alpha-y)^3]   ;

Fq1 = -inv(M)*[xi_h ; xi_alpha; 0];    % X^3
Fq2 =-inv(M)* [epsilon*xi ; -1*epsilon*xi*lambda ; -1*xi] ; % () 

 % 

 A = [zeros(3) , eye(3) ;
      -inv(M)*K , -inv(M) * C];
 


%
n=3;
dt = 0.1;
T_end = 100;
to=[0:dt:T_end];

x_o =  zeros(1,3);
x_o(1) = 1;
xdot_o =   zeros(1,3);
% xdot_o(1)=xx;
xo=[x_o xdot_o]; % Initial condition matrix for numerical simulation

%% Initial condition 
% x_o=[0.1*ones(n,1)]'; xdot_o=[0.1*ones(n,1)]'; % Initial displacements and velocities
E_o=.5*xdot_o*M*xdot_o'+.5*x_o*K*x_o'; % Initial energy

clear t;
% Numerical solution
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,x]=ode45(@odefun_66,to,xo,options,A,Fq1,Fq2); % x is the solution vector

% figure()
% plot(t,x(:,1),t,x(:,3),'LineWidth',2)
% hold on 
% grid on 
% % grid minor
% xlabel("Time (s)")
% ylabel("Response")
% % title("\gamma = \omega^2 / \omega_\alpha^2")
% % end
% 
% L = num2str(list' );
% legend(L)

%% Energy measures
YY=x(1:1+T_end/dt,1:n+1); % Responce corresponding to the transpose of x in the equations of motions
YY1=x(1:1+T_end/dt,n+1:2*n);%Responce corresponding to the transpose of x_dot in the equations of motions
NN=1+T_end/dt; % Number of time steps in the integration
Knl= [0;0;0;Fq1];
for i=1:NN
 Ed(i)=YY1(i:i,1:n)*C*(YY1(i:i,1:n))'; % Instantaneous energy dissipation by damping
 Eke(i)=.5*YY1(i:i,1:n)*M*(YY1(i:i,1:n))';% Instantaneous kinetic energy 
%  Epe_nL(i)=.25*Knl*(x(i,1)-x(i,2))^3; % Instantaneous nonlinear potential energy
%   Epe_nL(i)=.25*[0 ;0; 0; Fq1  ]'*x(i,:)'.^3+sum([0 ;0 ;0 ;Fq2].* (x(i,1)-x(i,2)-x(i,3))^3);
  Epe_nL(i)=.25* x(i,1:3)*(Fq1.*x(i,1:3)'.^3+Fq2.* (x(i,1)-x(i,2)-x(i,3))^3);
  Epe_L(i)=.5*YY(i:i,1:n)*K*(YY(i:i,1:n))'; % Instantaneous linear potential energy 
  Epe(i)=Epe_nL(i)+Epe_L(i);% Instantaneous total linear & nonlinear potential energy 
  Ef(i)=0;%Ef(i)=sum(YY1(i:i,1:n)*Ff*cos(w*t(i))); % Instantaneous forcing input energy 
end
%%
Elinear=cumtrapz(t,Ed); % Accumulated energy dissipated by linear damping 
Einput=cumtrapz(t,Ef);   % Accumulated energy input by forcing
figure()
subplot(2,2,1),plot(t,Eke+Epe,'k','LineWidth',2.5);hold on
xlabel('Time (s)','fontsize',16);
ylabel('E_p+E_k','fontsize',16);
hold on;grid on
set(gca,'fontsize',16)

subplot(2,2,2),plot(t,Elinear,'k','LineWidth',2.5);hold on
xlabel('Time (s)','fontsize',16);
ylabel('E_d_a_m_p_i_n_g','fontsize',16);
hold on;grid on
set(gca,'fontsize',16)

subplot(2,2,3),plot(t,Einput,'k','LineWidth',2.5);hold on
xlabel('Time (s)','fontsize',16);
ylabel('E_i_n_p_u_t','fontsize',16);
hold on;grid on
set(gca,'fontsize',16)

subplot(2,2,4),plot(t,Eke+Epe+Elinear,'k','LineWidth',2.5);hold on
subplot(2,2,4),plot(t,Einput+E_o,'b:','LineWidth',2.5);hold on
xlabel('Time (s)','fontsize',16);
ylabel('Energy (m)','fontsize',16);
hleg1=legend('E_k+E_p+E_d_a_m_p_i_n_g','E_i_n_p_u_t by forcing+E_o','FontSize',12);
set(hleg1, 'FontSize', 16);
set(gca,'fontsize',16)
hold on;grid on


