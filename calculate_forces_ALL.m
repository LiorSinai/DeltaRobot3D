function wrenches=calculate_forces_ALL(lambda,R_A1,R_A2,R_A3,R_L1,R_L2,R_L3,thetaA_t,L)
% This code is written for the controlled system
% Therefore the Jacobian used is Jacobian_d and the indices specfied
% below are for Qd, Phi_without_Driving and Jacobian_d
% If running a kinematically driven system, then run the script
% project_week_5.forces.m instead

N=size(thetaA_t,2); % number of valid points

%% Intialise the wrenches structure
% constraint forces for actuator 1
wrenches(1).point='A1'; % point
wrenches(1).body ='U1'; % reference body centroid. 
% Forces at the same point but viewed from different bodies should have the
% same magnitude but opposite directions
wrenches(1).indQd=[4:6,NaN]; % indices of the reference centroid in Q
wrenches(1).indPhi=[1:3,NaN]; % indices of the relevant constraint equations in Phi
r   =zeros(3,N);
for timeStep = 1:N
    R=R_A1*rot3D_Rodrigues([0;1;0],thetaA_t(1,timeStep));
    r(:,timeStep)=R*[0; 0 ;+L(1)/2]; 
end
wrenches(1).r=r; % vector from U1 to A1 in the inertial frame

wrenches(2).point='B1'; % point
wrenches(2).body ='U1'; % reference body centroid. 
wrenches(2).indQd=[4:6,NaN];
wrenches(2).indPhi=[4:6,NaN];
r=zeros(3,N); 
for timeStep = 1:N
    %R=eval(subs(R_U1,thetaA1,thetaA_t(1,timeStep)));
    R=R_A1*rot3D_Rodrigues([0;1;0],thetaA_t(1,timeStep));%=R_U1
    r(:,timeStep)=R*[0; 0 ;-L(1)/2]; % vector from U1 to B1  
end
wrenches(2).r=r;

% constraint forces for actuator 2
wrenches(3).point='A2'; % point
wrenches(3).body ='U2'; % reference body centroid. 
wrenches(3).indQd=[13:15,NaN];
wrenches(3).indPhi=[10:12,NaN];
r=zeros(3,N);
for timeStep = 1:N
    %R=eval(subs(R_U2,thetaA2,thetaA_t(2,timeStep))); %slow
    R=R_A2*rot3D_Rodrigues([0;1;0],thetaA_t(2,timeStep));%=R_U2
    r(:,timeStep)=R*[0; 0 ;+L(1)/2]; % vector from U2 to A2
end
wrenches(3).r=r;

wrenches(4).point='B2'; % point
wrenches(4).body ='U2'; % reference body centroid
wrenches(4).indQd=[13:15,NaN];
wrenches(4).indPhi=[13:15,NaN];
r=zeros(3,N); 
for timeStep = 1:N
    %R=eval(subs(R_U2,thetaA2,thetaA_t(2,timeStep))); %slow
    R=R_A2*rot3D_Rodrigues([0;1;0],thetaA_t(2,timeStep));
    r(:,timeStep)=R*[0; 0 ;-L(1)/2]; % vector from U2 to B2 
end
wrenches(4).r=r;

% constraint forces for actuator 3
wrenches(5).point='A3'; % point
wrenches(5).body ='U3'; % reference body centroid
wrenches(5).indQd=[22:24,NaN];
wrenches(5).indPhi=[19:21,NaN];
r=zeros(3,N);
for timeStep = 1:N
    %R=eval(subs(R_U3,thetaA3,thetaA_t(3,timeStep))); %slow
    R=R_A3*rot3D_Rodrigues([0;1;0],thetaA_t(3,timeStep));
    r(:,timeStep)=R*[0; 0 ;+L(1)/2]; % vector from U3 to A3
end
wrenches(5).r=r;

wrenches(6).point='B3'; % point
wrenches(6).body ='U3'; % reference body centroid
wrenches(6).indQd=[22:24,NaN];
wrenches(6).indPhi=[22:24,NaN];
r=zeros(3,N); 
for timeStep = 1:N
    %R=eval(subs(R_U3,thetaA3,thetaA_t(3,timeStep))); %slow
    R=R_A3*rot3D_Rodrigues([0;1;0],thetaA_t(3,timeStep));
    r(:,timeStep)=R*[0; 0 ;-L(1)/2]; % vector from U3 to B3  
end
wrenches(6).r=r;

%% Reaction forces
tic
fprintf('Calculating reaction forces...\n')
for currentPoint=1:length(wrenches)
    fprintf('%d point(s) complete of %d points\n',currentPoint-1,length(wrenches))
    [wrenches(currentPoint).forces,wrenches(currentPoint).moments] = ...
        calculate_forces_dependent(R_L1,R_L2,R_L3,thetaA_t,L,lambda,...
        wrenches(currentPoint).indQd,wrenches(currentPoint).indPhi,...
        wrenches(currentPoint).r);
end
fprintf('%d point(s) complete of %d points\n',currentPoint,length(wrenches))
toc