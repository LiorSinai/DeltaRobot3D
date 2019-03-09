%% Test simulations
% Lior Sinai, 2019-02-23
% run project_week_5.m and controller_Delta3D_inverse.slx first.

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
Qacc_sim(ind_d,:)=calculate_acceleration_dependent(...
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
% map to kinematically driven variable names
thetaA_t=thetaA_sim;
Qacc=Qacc_sim;
time_range=t_sim;
Qvel=Qvel_sim;
wrenches(1).momentsC=[0;1;0].*tau_sim(1,:);
wrenches(3).momentsC=[0;1;0].*tau_sim(2,:);
wrenches(5).momentsC=[0;1;0].*tau_sim(3,:);

% now can run the old test_force.m
run test_forces.m
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
