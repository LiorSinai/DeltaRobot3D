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

%% plot desired trajectories
figure;
% show full trajectory
plot_boundaries(L,thetaU1_z_0,thetaU2_z_0,thetaU3_z_0) 
th=linspace(0,2*pi,50);
pRef=[cos(3*th);sin(3 *th);-1.8+0.5*th];
plot3(pRef(1,:), pRef(2,:), pRef(3,:))
hold on
plot3(pRef(1,1), pRef(2,1), pRef(3,1),'x')

%% Plot simulations
plot_Delta3D( Q_sim,L,R_A1,R_A2,R_A3,R_L1_sim,R_L2_sim,R_L3_sim,...
    thetaA_sim,t_sim,0.2,false )

