%% Test simulations
% Lior Sinai, 2019-02-23

SIZE_QD=33;
SIZE_QI=9;
SIZE_Q=42;
%% extract data
% initialise
N=length(Qi_SIM.Time);
Q_sim=zeros(SIZE_Q,N);
Qvel_sim=Q_sim;
Qacc_sim=Q_sim;

% independent and dependent co-ordinate indices
ind_i=[7:9,19:21,31:33]';
ind_d=(1:42)';
ind_d(ind_i)=[]; % remove independent co-ordinates

% set data
Q_sim(ind_i,:)= Qi_SIM.Data';
Q_sim(ind_d,:)= Qd_SIM.Data';
Qvel_sim(ind_i,:)=Qi_vel_SIM.Data';
Qvel_sim(ind_d,:)=Qd_vel_SIM.Data';
Qacc_sim(ind_i,:)=Qi_acc_SIM.Data';

t_sim = Qi_SIM.Time;
thetaA_sim = thetaA_SIM.Data';
%tau_sim=tau_fb_SIM.Data'; 
tau_sim=tau_a_SIM.Data';
%tau_sim=zeros(size(tau_a_SIM.Data'));

R_L1_sim=zeros(3,3,N);
R_L2_sim=R_L1_sim;
R_L3_sim=R_L1_sim;
for timeStep=1:N
    R_L1_sim(:,:,timeStep)=R_L1_SIM.Data(:,:,timeStep);
    R_L2_sim(:,:,timeStep)=R_L2_SIM.Data(:,:,timeStep);
    R_L3_sim(:,:,timeStep)=R_L3_SIM.Data(:,:,timeStep);
end

%% Check actuator angles
%Simple test to check that the actuator angles are being calculated correctly
% vecnorm(Q_sim(7:8,:)-  [0;thetaU1_y_0])-abs(thetaA_sim(1,:)-thetaA_0(1))
% vecnorm(Q_sim(19:20,:)-[0;thetaU2_y_0])-abs(thetaA_sim(2,:)-thetaA_0(2))
% vecnorm(Q_sim(31:32,:)-[0;thetaU3_y_0])-abs(thetaA_sim(3,:)-thetaA_0(3))

%% Check deviations of the rotation matrices
detR1=zeros(size(t_sim,2),1);
detR2=detR1;
detR3=detR1;
% The determinant should be 1 at all time steps
for timeStep=1:N
    detR1(timeStep)=det(R_L1_sim(:,:,timeStep)); 
    detR2(timeStep)=det(R_L2_sim(:,:,timeStep));  
    detR3(timeStep)=det(R_L3_sim(:,:,timeStep));  
end
figure;
plot(1:N,detR1-1,'-x',1:N,detR2-1,'-*',1:N,detR3-1,'-o')
legend('det(R_{L1})','det(R_{L2})','det(R_{L3})')
xlabel('time steps')
title('Error of rotation matrices')

%% calculate the forces
% depedent accelerations are not calculated during simulation, so they have
% to be estimated
Qacc_sim(ind_d,:)=caculate_acceleration_dependent(...
    Qacc_sim(ind_i,:),Qvel_sim,...
    R_L1_sim,R_L2_sim,R_L3_sim,thetaA_sim,...
    L);
[lambda,Q_Ad]=calculate_lagrange_dependent(...
    Qacc_sim(ind_i,:),Qacc_sim(ind_d,:),Qvel_sim,...
    R_A1,R_A2,R_A3,...
    R_L1_sim,R_L2_sim,R_L3_sim,thetaA_sim,...
    L,I_C,m_u,m_l,m_p,g);
wrenches=calculate_forces_ALL(lambda,...
    R_A1,R_A2,R_A3,...
    R_L1_sim,R_L2_sim,R_L3_sim,thetaA_sim,...
    L);

%% Check these forces
% note the moments are different test_forces.m
errorForcesResidue=zeros(9,N);
errorMomentsResidue=zeros(9,N);
clear M1 M2 M3 M4 M5
for timeStep = 1:N
    % map values
    forcesA1=wrenches(1).forces(:,timeStep);
    forcesB1=wrenches(2).forces(:,timeStep);
    forcesA2=wrenches(3).forces(:,timeStep);
    forcesB2=wrenches(4).forces(:,timeStep);
    forcesA3=wrenches(5).forces(:,timeStep);
    forcesB3=wrenches(6).forces(:,timeStep);
    
    % force balance in the inertial frame 0
    % m_u*rddot =FA+FB+[0;0;-m_u*g]
    errorForcesResidue(1:3,timeStep)=m_u*Qacc_sim(4:6,timeStep)+...
                               -forcesA1-forcesB1-[0;0;-m_u*g];
   errorForcesResidue(4:6,timeStep)=m_u*Qacc_sim(16:18,timeStep)+...
                               -forcesA2-forcesB2-[0;0;-m_u*g];
   errorForcesResidue(7:9,timeStep)=m_u*Qacc_sim(28:30,timeStep)+...
                               -forcesA3-forcesB3-[0;0;-m_u*g];
   
   % moment balance in the body fixed frame C
   % I_c*alphaC+tilde(omegaC)*I_c*omegaC = tilde(r^U_A)*F^A_U+tilde(r^U_B)*F^A_B
   %                                       +[0;M_A;0]
   %R=eval(subs(R_U1,thetaA1,thetaA_t(1,timeStep)));
   R=R_A1*rot3D_Rodrigues([0;1;0],thetaA_sim(1,timeStep));%=R_U1 faster
   alphaC=R'*Qacc_sim(7:9,timeStep);
   omegaC=R'*Qvel_sim(7:9,timeStep); 
   momentsC=[0;tau_sim(1,timeStep);0]; % independent moment
   errorMomentsResidue(1:3,timeStep)=I_C(1:3).*alphaC+...
                                +tilde(omegaC)*diag(I_C(1:3))*omegaC+...
                                -tilde([0; 0 ;+L(1)/2])*R'*forcesA1+...
                                -tilde([0; 0 ;-L(1)/2])*R'*forcesB1+...
                                -momentsC;
   % components of the error
   M1(:,timeStep)=I_C(1:3).*alphaC;
   M2(:,timeStep)=tilde(omegaC)*diag(I_C(1:3))*omegaC;
   M3(:,timeStep)=tilde([0; 0 ;+L(1)/2])*R'*forcesA1;
   M4(:,timeStep)=tilde([0; 0 ;-L(1)/2])*R'*forcesB1;
   M5(:,timeStep)=momentsC;
   textLegend={'$I\alpha$','','',...
               '$\tilde{\omega}I\omega$','','',...
               '$r_A\times R^T F_{A,1}$','','',...
               '$r_B\times R^T F_{B,1}$','','',...
               '$M_{a}$','',''};
   
   %R=eval(subs(R_U2,thetaA2,thetaA_t(2,timeStep)));  % slow
   R=R_A2*rot3D_Rodrigues([0;1;0],thetaA_sim(2,timeStep));%=R_U2
   omegaC=R'*Qvel_sim(19:21,timeStep); 
   alphaC=R'*Qacc_sim(19:21,timeStep);
   momentsC=[0;tau_sim(2,timeStep);0];
   errorMomentsResidue(4:6,timeStep)=I_C(1:3).*alphaC+...
                                +tilde(omegaC)*diag(I_C(1:3))*omegaC+...
                                -tilde([0; 0 ;+L(1)/2])*R'*forcesA2+...
                                -tilde([0; 0 ;-L(1)/2])*R'*forcesB2+...
                                -momentsC;

   %R=eval(subs(R_U3,thetaA3,thetaA_t(3,timeStep)));
   R=R_A3*rot3D_Rodrigues([0;1;0],thetaA_sim(3,timeStep)); %=R_U3
   omegaC=R'*Qvel_sim(31:33,timeStep); 
   alphaC=R'*Qacc_sim(31:33,timeStep);
   momentsC=[0;tau_sim(3,timeStep);0];
   errorMomentsResidue(7:9,timeStep)=I_C(1:3).*alphaC+...
                                +tilde(omegaC)*diag(I_C(1:3))*omegaC+...
                                -tilde([0; 0 ;+L(1)/2])*R'*forcesA3+...
                                -tilde([0; 0 ;-L(1)/2])*R'*forcesB3+...
                                -momentsC;
                            
end

%% plot components
% Note: for twist, there will be an applied moment about both the body frame
% x-axis and y-axis for arms 2 and 3.
% This does not seem to be an error, as the asymmetry of the sitatuion
% tends to rotate the delta robot about the x axis, and therefore
% an applied moment has to counter-act this
figure;plot(1:N,M1,'-xb',1:N,M2,'-xr',1:N,M3,'-xg',...
            1:N,M4,'-xm',1:N,M5,'-xk');
title('Components of the forces')
grid;legend(textLegend,'interpreter','latex');

%% plot error residues
figure; 
plot(1:N,errorForcesResidue,'-x'); legend; grid
title('Error residues for forces')
figure; 
plot(1:N,errorMomentsResidue,'-x'); legend; grid
title('Error residues for moments')

%save(sprintf('TestRun_%s',timestamp)); % save all variables

%% Test inverse kinematics
% qa_t=zeros(size(thetaA_sim));
% qaVel_t=qa_t;
% qaAcc_t=qa_t;
% for testIndex=1:N
%     qa0=thetaA_sim(:,testIndex);
%     [qa, qaVel, qaAcc]=calculate_IK(L,qa0,Q_sim(indP,testIndex),...
%         Qvel_sim(indP,testIndex),Qacc_sim(indP,testIndex),...
%         tolerance);
%     qa_t(:,testIndex)=qa;
%     qaVel_t(:,testIndex)=qaVel;
%     qaAcc_t(:,testIndex)=qaAcc;
% end
% 
% figure;
% plot(1:N,qa_t-thetaA_sim,'x-')
% title('Error in \theta_A')
% 
% % figure;
% % plot(1:N,qaVel_t-omega,'x-')
% % title('Error in \omega_A')
% 
% figure;
% plot(1:N,qaAcc_t,'x-')
% title('Error in \alpha_A')
