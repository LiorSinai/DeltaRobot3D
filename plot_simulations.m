% Plot simulation results
N=length(Qi_SIM.Time);
Q_sim=zeros(42,N);

% independent and dependent co-ordinate indices
ind_i=[7:9,19:21,31:33]';
ind_d=(1:42)';
ind_d(ind_i)=[]; % remove independent co-ordinates

Q_sim(ind_i,:)= Qi_SIM.Data';
Q_sim(ind_d,:)= Qd_SIM.Data';
t_sim = Qi_SIM.Time;
thetaA_sim = thetaA_SIM.Data';

R_L1_sim=zeros(3,3,N);
R_L2_sim=R_L1_sim;
R_L3_sim=R_L1_sim;
for timeStep=1:N
    R_L1_sim(:,:,timeStep)=R_L1_SIM.Data(:,:,timeStep);
    R_L2_sim(:,:,timeStep)=R_L2_SIM.Data(:,:,timeStep);
    R_L3_sim(:,:,timeStep)=R_L3_SIM.Data(:,:,timeStep);
end
plot_Delta3D( Q_sim,L,R_A1,R_A2,R_A3,R_sym,thetaA_sym,R_L1_sim,R_L2_sim,R_L3_sim,thetaA_sim,t_sim,0.2 )

%% check simulations
% vecnorm(Q_sim(7:8,:)-  [0;thetaU1_y_0])-abs(thetaA_sim(1,:)-thetaA_0(1))
% vecnorm(Q_sim(19:20,:)-[0;thetaU2_y_0])-abs(thetaA_sim(2,:)-thetaA_0(2))
% vecnorm(Q_sim(31:32,:)-[0;thetaU3_y_0])-abs(thetaA_sim(3,:)-thetaA_0(3))

