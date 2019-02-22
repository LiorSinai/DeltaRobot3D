% 3 February 2019
% Test Embedded Formulation
% Should run prokect_week_5/m first to get:
%   all symbols
%   all rotation matrices
%   Q_AC

% independent co-ordinates  
% Qi = [thetaU1_x;thetaU1_y;thetaU1_z;
%       thetaU2_x;thetaU2_y;thetaU2_z;
%       thetaU3_x;thetaU3_y;thetaU3_z;]; 
ind_i=[7:9,19:21,31:33]';
% dependent co-ordinates
% Qd=[M_x;  M_y;  M_z; 
%     % arm 1
%    U1_x; U1_y; U1_z; 
%    L1_x; L1_y; L1_z; thetaL1_x; thetaL1_y; thetaL1_z;
%     % arm 2
%    U2_x; U2_y; U2_z; 
%    L2_x; L2_y; L2_z; thetaL2_x; thetaL2_y; thetaL2_z;
%     % arm 3
%    U3_x; U3_y; U3_z; 
%    L3_x; L3_y; L3_z; thetaL3_x; thetaL3_y; thetaL3_z;
%    P_x;  P_y; P_z
%   ];
ind_d=(1:42)';
ind_d(ind_i)=[]; % remove independent co-ordinates


Jacobian_without_driving=[Jacobian(1:9,:) ; Jacobian(13:21,:) ; Jacobian(25:33,:); Jacobian(37:42,:)];
Gamma_without_driving=[Gamma(1:9,:) ; Gamma(13:21,:) ; Gamma(25:33,:); Gamma(37:42,:)];


% Calculate moment inertias in the inertia frame
I_U1_sym=R_U1*diag([I_upper_xx,I_upper_yy, I_upper_zz])*R_U1.';
I_U2_sym=R_U2*diag([I_upper_xx,I_upper_yy, I_upper_zz])*R_U2.';
I_U3_sym=R_U3*diag([I_upper_xx,I_upper_yy, I_upper_zz])*R_U3.';

I_L1=R_L1*diag([I_lower_xx,I_lower_yy, I_lower_zz])*R_L1.';
I_L2=R_L2*diag([I_lower_xx,I_lower_yy, I_lower_zz])*R_L2.';
I_L3=R_L3*diag([I_lower_xx,I_lower_yy, I_lower_zz])*R_L3.';

% initialise acceleration matrix
N=size(Q,2);
Qacc_est=zeros(9,N);
GammEst=zeros(33,N);
for timeStep=1:N % test index
    fprintf('% d ... %.1f%% complete of testing\n',timeStep,100*timeStep/length(time_range))
    % Calculate moment inertias in the inertia frame

    I_U1_t=eval(subs(I_U1_sym,thetaA_sym(1),thetaA_t(1,timeStep)));
    I_U2_t=eval(subs(I_U2_sym,thetaA_sym(2),thetaA_t(2,timeStep)));
    I_U3_t=eval(subs(I_U3_sym,thetaA_sym(3),thetaA_t(3,timeStep)));

    I_L1_t=eval(subs(I_L1,R_L1,R1(:,:,timeStep)));
    I_L2_t=eval(subs(I_L2,R_L2,R2(:,:,timeStep)));
    I_L3_t=eval(subs(I_L3,R_L3,R3(:,:,timeStep)));

    Mii = [I_U1_t zeros(3,6);
          zeros(3,3) I_U2_t  zeros(3,3);
          zeros(3,6) I_U3_t];

    Mdd = diag([0 0 0 ...
               m_u*[1 1 1]              ...
               m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
               m_u*[1 1 1]              ...
               m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
               m_u*[1 1 1]              ...
               m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
               m_p*[1 1 1] ]); 
    % Correct Mdd
    Mdd(10:12,10:12)=I_L1_t;
    Mdd(19:21,19:21)=I_L2_t;
    Mdd(28:30,28:30)=I_L3_t;

    Mdi = zeros(33,9);
    Mid = zeros(9,33);

    Q_Ad=Q_AC(ind_d,timeStep);
    Q_Ai=Q_AC(ind_i,timeStep);
    Q_Ai(1:3)=Q_Ai(1:3)+wrenches(1).moments(:,timeStep);
    Q_Ai(4:6)=Q_Ai(4:6)+wrenches(3).moments(:,timeStep);
    Q_Ai(7:9)=Q_Ai(7:9)+wrenches(5).moments(:,timeStep);

    Jcalc=subs(Jacobian_without_driving,[R_L1,R_L2,R_L3],[R1(:,:,timeStep),R2(:,:,timeStep),R3(:,:,timeStep)]);
    Jcalc=subs(Jcalc,thetaA_sym,thetaA_t(:,timeStep));
    Jcalc=eval(subs(Jcalc,Q_sym,Q(:,timeStep)));

    G=subs(Gamma_without_driving,t_sym,time_range(timeStep));
    G=subs(G,[R_L1,R_L2,R_L3],[R1(:,:,timeStep),R2(:,:,timeStep),R3(:,:,timeStep)]);
    G=subs(G,thetaA_sym,thetaA_t(:,timeStep));
    G=subs(G,Qvel_sym,Qvel(:,timeStep));
    G=eval(subs(G,Q_sym,Q(:,timeStep)));
    
    %Jcalc*Qddot(:,step)-G

    Ji=Jcalc(:,ind_i);
    Jd=Jcalc(:,ind_d);

    J =Jd\Ji;
    J2=Jd\G;

    Mhat = Mii - Mid*J - J.'*(Mdi - Mdd*J);
    Qhat =Q_Ai-Mid*J2- (J.')*(Q_Ad - Mdd*J2);
    
    GammEst(:,timeStep)=G;
    
    Qacc_est(:,timeStep)=Mhat\Qhat;
    %Qddot(ind_i,step); % these should be equal
end

%% Plot results
figure;
plot(1:N,Qacc(ind_i,:),'-',1:N,Qacc_est,'-x')
title('Prescribed vs estimated accelerations (symbolic)')

%% Test simulink function
addpath('Numeric Functions')

% initialise acceleration matrices
N=size(Q,2);
Qacc_est2=zeros(9,N);
Qvel_est2 =zeros(42,N);
GammaEst2 =zeros(33,N);
QHAT2=zeros(9,N);
for timeStep=1:N % test index
    fprintf('% d ... %.1f%% complete of testing\n',timeStep,100*timeStep/length(time_range))
    
    qi=Q(ind_i,timeStep);
    qi_vel=Qvel(ind_i,timeStep);
    
    [Mhat, Qhat,qd_vel_est,G_est2,qdot] = ...
    computeMhatQhat(qi,qi_vel,R1(:,:,timeStep),R2(:,:,timeStep),R3(:,:,timeStep),...
        m_u, m_p, m_l, L, g,...
        I_upper_xx, I_upper_yy, I_upper_zz,...
        I_lower_xx, I_lower_yy, I_lower_zz,...
        thetaA_0);
    
    GammaEst2(:,timeStep)=G_est2;
    
    
    Qhat(1:3)=Qhat(1:3)+wrenches(1).moments(:,timeStep);
    Qhat(4:6)=Qhat(4:6)+wrenches(3).moments(:,timeStep);
    Qhat(7:9)=Qhat(7:9)+wrenches(5).moments(:,timeStep);
    
    QHAT2(:,timeStep)=Qhat;

    Qvel_est2(ind_i,timeStep)=qi_vel;
    Qvel_est2(ind_d,timeStep)=qd_vel_est;
    
    Qacc_est2(:,timeStep)=Mhat\Qhat;
    %Qddot(ind_i,timeStep); % these should be equal
end

%% plot results

figure;
plot(1:N,Qacc(ind_i,:),'-',1:N,Qacc_est2,'-x')
title('Prescribed vs estimated accelerations (numeric)')

figure;
plot(1:N,Qvel_est2(ind_d,:)-Qvel(ind_d,:))
title('Errors for dependent velocities')

figure;

plot(1:N,GammaEst2-GammEst)
% plot(1:N,Q_est(ind_d,:)-Q_AC(ind_d,:))
title('Comparison of \gamma for symbolic and numeric calculations')
