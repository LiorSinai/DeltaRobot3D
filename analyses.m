%% Analyse results
N=length(Qi_CTC.Time);

detR1=zeros(size(t_sim,2),1);
detR2=detR1;
detR3=detR1;
for timeStep=1:N
    detR1(timeStep)=det(R_L1_sim(:,:,timeStep));    
    detR2(timeStep)=det(R_L2_sim(:,:,timeStep));  
    detR3(timeStep)=det(R_L3_sim(:,:,timeStep));  
end
figure;
plot(1:N,detR1,'x',1:N,detR2,'*',1:N,detR3,'o')
legend('det(R_{L1})','det(R_{L2})','det(R_{L3})')
xlabel('time steps')
title('Growth in error of rotation matrices')

% check errors
% all accelerations at maximum/minimum velocities = 0
%% Get velocities, accelerations and tangential accelerations
U1v=velocities(4:6,:);
U1a=accelerations(4:6,:);
u=U1v./vecnorm(U1v); % unit vector in the tangential direction
U1a_tan=sum(u.*U1a); % dot product

omegaU1=velocities(7:9,:);
alphaU1=accelerations(7:9,:);
u=omegaU1./vecnorm(omegaU1); % unit vector of axis of rotation
alphaU1_tan=sum(u.*alphaU1);

L1v=velocities(10:12,:);
L1a=accelerations(10:12,:);
u=L1v./vecnorm(L1v); % unit vector in the tangential direction
L1a_tan=sum(u.*L1a); % dot product

omegaL1=velocities(13:15,:);
alphaL1=accelerations(13:15,:);
u=omegaL1./vecnorm(omegaL1); % unit vector of axis of rotation
alphaL1_tan=sum(u.*alphaL1);

% second center of mass     
U2v=velocities(16:18,:);
U2a=accelerations(16:18,:);
u=U2v./vecnorm(U2v); % unit vector in the tangential direction
U2a_tan=sum(u.*U2a); % dot product

omegaU2=velocities(19:21,:);
alphaU2=accelerations(19:21,:);
u=omegaU2./vecnorm(omegaU2); % unit vector of axis of rotation
alphaU2_tan=sum(u.*alphaU2);

L2v=velocities(22:24,:);
L2a=accelerations(22:24,:);
u=L2v./vecnorm(L2v); % unit vector in the tangential direction
L2a_tan=sum(u.*L2a); % dot product

omegaL2=velocities(25:27,:);
alphaL2=accelerations(25:27,:);
u=omegaL2./vecnorm(omegaL2); % unit vector of axis of rotation
alphaL2_tan=sum(u.*alphaL2);

% third center of mass
U3v=velocities(28:30,:);
U3a=accelerations(28:30,:);
u=U3v./vecnorm(U3v); % unit vector in the tangential direction
U3a_tan=sum(u.*U3a); % dot product

omegaU3=velocities(31:33,:);
alphaU3=accelerations(31:33,:);
u=omegaU3./vecnorm(omegaU3); % unit vector of axis of rotation
alphaU3_tan=sum(u.*alphaU3);

L3v=velocities(34:36,:);
L3a=accelerations(34:36,:);
u=L3v./vecnorm(L3v); % unit vector in the tangential direction
L3a_tan=sum(u.*L3a); % dot product

omegaL3=velocities(37:39,:);
alphaL3=accelerations(37:39,:);
u=omegaL3./vecnorm(omegaL3); % unit vector of axis of rotation
alphaL3_tan=sum(u.*alphaL3);

%end-effector

Pv = velocities(40:42,:);
Pa = accelerations(40:42,:);
u=Pv./vecnorm(Pv); % unit vector in the tangential direction
Pa_tan=sum(u.*Pa); % dot product

%% Check accelerations at mins/maximums
% At max(velocity),acceleration=0
errors=zeros(12,1);
[~,ind]=max(vecnorm(U1v));
errors(1)=U1a_tan(ind)./max(abs(U1a_tan)); 
[~,ind]=max(vecnorm(omegaU1));
errors(2)=alphaU1_tan(ind)./max(abs(alphaU1_tan)); 

[~,ind]=max(vecnorm(L1v));
errors(3)=L1a_tan(ind)./max(abs(L1a_tan));
[~,ind]=max(vecnorm(omegaL1));
errors(4)=alphaL1_tan(ind)./max(abs(alphaL1_tan)); 

[~,ind]=max(vecnorm(U2v));
errors(5)=U2a_tan(ind)./max(abs(U2a_tan));
[~,ind]=max(vecnorm(omegaU2));
errors(6)=alphaU2_tan(ind)./max(abs(alphaU2_tan));

[~,ind]=max(vecnorm(L2v));
errors(7)=L2a_tan(ind)./max(abs(L2a_tan));
[~,ind]=max(vecnorm(omegaL2));
errors(8)=alphaL2_tan(ind)./max(abs(alphaL2_tan));

[~,ind]=max(vecnorm(U3v));
errors(9)=U3a_tan(ind)./max(abs(U3a_tan));
[~,ind]=max(vecnorm(omegaU3));
errors(10)=alphaU3_tan(ind)./max(abs(alphaU3_tan));

[~,ind]=max(vecnorm(L3v));
errors(11)=L3a_tan(ind)./max(abs(L3a_tan));
[~,ind]=max(vecnorm(omegaL3));
errors(12)=alphaL3_tan(ind)./max(abs(alphaL3_tan)); 

[~,ind]=max(vecnorm(Pv));
errors(13)=Pa_tan(ind)./max(abs(Pa_tan));  

figure;
plot(1:length(errors),errors,'x')

