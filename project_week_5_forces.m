%% caclulate forces
% continuation from project_week_5.m
% run project_week_5.m first with calculate=true

% Calculates forces and moments base on the Lagrange multiplier method 
% tic
% [lambda, Q_AC]=calculate_lagrange(Jacobian,Q_sym,R_sym,thetaA_sym,Q,Qvel,Qacc,R1,R2,R3,thetaA_t,time_range,M_C,I_C,Q_A);
% toc
%tic
[lambda, Q_AC] = calculate_lagrange_fast(Q,Qvel,Qacc,R1,R2,R3,thetaA_t,time_range,L,M_C,I_C,Q_A);
%toc

Nvalid=size(thetaA_t,2); % number of valid points
time_valid=time_range(1:Nvalid); % if a singularity was hit, the valid_time
                                 % will be shorter than the time_range 

%% Intialise the wrenches structure
% constraint forces for actuator 1
wrenches(1).point='A1'; % point
wrenches(1).body ='U1'; % reference body centroid. 
% Forces at the same point but viewed from different bodies should have the
% same magnitude but opposite directions
wrenches(1).indQ=4:9; % indices of the reference centroid in Q
wrenches(1).indPhi=[1:3,10:12]; % indices of the relevant constraint equations in Phi
r   =zeros(3,Nvalid);
for timeStep = 1:Nvalid
    R=R_A1*rot3D_Rodrigues([0;1;0],thetaA_t(1,timeStep));
    r(:,timeStep)=R*[0; 0 ;+L_u/2]; 
end
wrenches(1).r=r; % vector from U1 to A1 in the inertial frame

wrenches(2).point='B1'; % point
wrenches(2).body ='U1'; % reference body centroid. 
wrenches(2).indQ=4:9;
wrenches(2).indPhi=[4:6,10:12];
r=zeros(3,Nvalid); 
for timeStep = 1:Nvalid
    %R=eval(subs(R_U1,thetaA1,thetaA_t(1,timeStep)));
    R=R_A1*rot3D_Rodrigues([0;1;0],thetaA_t(1,timeStep));%=R_U1
    r(:,timeStep)=R*[0; 0 ;-L_u/2]; % vector from U1 to B1  
end
wrenches(2).r=r;

% constraint forces for actuator 2
wrenches(3).point='A2'; % point
wrenches(3).body ='U2'; % reference body centroid. 
wrenches(3).indQ=16:21;
wrenches(3).indPhi=[13:15,22:24];
r=zeros(3,Nvalid);
for timeStep = 1:Nvalid
    %R=eval(subs(R_U2,thetaA2,thetaA_t(2,timeStep))); %slow
    R=R_A2*rot3D_Rodrigues([0;1;0],thetaA_t(2,timeStep));%=R_U2
    r(:,timeStep)=R*[0; 0 ;+L_u/2]; % vector from U2 to A2
end
wrenches(3).r=r;

wrenches(4).point='B2'; % point
wrenches(4).body ='U2'; % reference body centroid
wrenches(4).indQ=16:21;
wrenches(4).indPhi=[16:18,22:24];
r=zeros(3,Nvalid); 
for timeStep = 1:Nvalid
    %R=eval(subs(R_U2,thetaA2,thetaA_t(2,timeStep))); %slow
    R=R_A2*rot3D_Rodrigues([0;1;0],thetaA_t(2,timeStep));
    r(:,timeStep)=R*[0; 0 ;-L_u/2]; % vector from U2 to B2 
end
wrenches(4).r=r;

% constraint forces for actuator 3
wrenches(5).point='A3'; % point
wrenches(5).body ='U3'; % reference body centroid
wrenches(5).indQ=28:33;
wrenches(5).indPhi=[25:27,34:36];
r=zeros(3,Nvalid);
for timeStep = 1:Nvalid
    %R=eval(subs(R_U3,thetaA3,thetaA_t(3,timeStep))); %slow
    R=R_A3*rot3D_Rodrigues([0;1;0],thetaA_t(3,timeStep));
    r(:,timeStep)=R*[0; 0 ;+L_u/2]; % vector from U3 to A3
end
wrenches(5).r=r;

wrenches(6).point='B3'; % point
wrenches(6).body ='U3'; % reference body centroid
wrenches(6).indQ=28:33;
wrenches(6).indPhi=[28:30,34:36];
r=zeros(3,Nvalid); 
for timeStep = 1:Nvalid
    %R=eval(subs(R_U3,thetaA3,thetaA_t(3,timeStep))); %slow
    R=R_A3*rot3D_Rodrigues([0;1;0],thetaA_t(3,timeStep));
    r(:,timeStep)=R*[0; 0 ;-L_u/2]; % vector from U3 to B3  
end
wrenches(6).r=r;

%% Reaction forces
tic
fprintf('Calculating reaction forces...\n')
for currentPoint=1:length(wrenches)
    fprintf('%d point(s) complete of %d points\n',currentPoint-1,length(wrenches))
%     [wrenches(currentPoint).forces,wrenches(currentPoint).moments] = ...
%         calculate_forces(Jacobian,Q_sym,R_sym,thetaA_sym,...
%         Q,R1,R2,R3,thetaA_t,time_range,lambda,...
%         wrenches(currentPoint).indQ,wrenches(currentPoint).indPhi,...
%         wrenches(currentPoint).r);
    [wrenches(currentPoint).forces,wrenches(currentPoint).moments] = ...
        calculate_forces_fast(R1,R2,R3,thetaA_t,time_range,L,lambda,...
        wrenches(currentPoint).indQ,wrenches(currentPoint).indPhi,...
        wrenches(currentPoint).r);
end
fprintf('%d point(s) complete of %d points\n',currentPoint,length(wrenches))
toc
% moments at B should be equal to moments at A

%% Transfer moments to body centre frames

% actuator 1
wrenches(1).momentsC=zeros(size(wrenches(1).moments));
% actuator 2
wrenches(3).momentsC=wrenches(1).momentsC;
% actuator 3
wrenches(5).momentsC=wrenches(1).momentsC;

for timeStep = 1:Nvalid
    %R=eval(subs(R_U1,thetaA1,thetaA_t(1,timeStep))); %slow
    R=R_A1*rot3D_Rodrigues([0;1;0],thetaA_t(1,timeStep));
    wrenches(1).momentsC(:,timeStep)=R'*wrenches(1).moments(:,timeStep);
    
    %R=eval(subs(R_U2,thetaA2,thetaA_t(2,timeStep))); %slow
    R=R_A2*rot3D_Rodrigues([0;1;0],thetaA_t(2,timeStep));
    wrenches(3).momentsC(:,timeStep)=R'*wrenches(3).moments(:,timeStep);
    
    %R=eval(subs(R_U3,thetaA3,thetaA_t(3,timeStep))); % slow
    R=R_A3*rot3D_Rodrigues([0;1;0],thetaA_t(3,timeStep));
    wrenches(5).momentsC(:,timeStep)=R'*wrenches(5).moments(:,timeStep);
end

%%
figure('Name','Moments arm 3');
hold on
colors={'k','b','r'};
for k=1:3
    plot(time_valid,wrenches(5).moments(k,:),'-','Color',colors{k})
end
hold on
for k=1:3
    plot(time_valid,wrenches(5).momentsC(k,:),'x','Color',colors{k})
end
title('Comparison of moments in the inertial and body fixed frame')
legend('M3_x','M3_y','M3_z','MC3_x','MC3_y','MC3_z')

% For PushDown, these should all align
figure('Name','Moments');
specs={'-','-.','--x'};
hold on
for k=1:3
    plot(time_valid,wrenches(1).momentsC(k,:),specs{k},'Color',colors{1})
end
for k=1:3
    plot(time_valid,wrenches(3).momentsC(k,:),specs{k},'Color',colors{3})
end
for k=1:3
    plot(time_valid,wrenches(5).momentsC(k,:),specs{k},'Color',colors{2})
end
title('Moments in the body fixed frames')
legend('MC1_x','MC1_y','MC1_z', ...
       'MC2_x','MC2_y','MC2_z', ... 
       'MC3_x','MC3_y','MC3_z' ...
       );

%save(sprintf('TestRun_%s',timestamp)); % save all variables
