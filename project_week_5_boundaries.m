% plot_inverse
% Calculates and plots the boundaries for a Delta Robot
clear

SIZE_Q=42;
timestamp=datestr(now,'yyyymmdd_HHMMSS');
addpath('Numeric Functions')

% input parameters
L = [0.8,1,0.4,0.5]; % L_upper, L_lower, L_effector, L_base
%L =[0.8;0.8; 0.6;0.6];
tolerance=0.0001;

%% independent and dependent co-ordinate indices
ind_i=[7:9,19:21,31:33]';
ind_d=(1:42)';
ind_d(ind_i)=[]; % remove independent co-ordinates

%% Set a good initial position
%%% Important for Simulink
%%%             z
%%%             |
%%% Start at 1__M__2_ _ _x          y
%%%          /     \U           3\ |
%%%          \     /L             \|_ _2 x
%%%           \_P_/               /
%%%                             1/
q0=ones(SIZE_Q,1);
%arm1
%arm1
%thetaU1_x_0=0; not free to choose. This will cause problems if this is not
%zero, because then R(1,2) and R(2,2) will be a function of theta_U1_y
%See the driving constraints
thetaU1_y_0=-pi/4;  
thetaU1_z_0=-2*pi/3; % angle in the x-y plane
% thetaL1_x_0=0;
thetaL1_y_0=pi/5;
thetaL1_z_0=thetaU1_z_0;
%arm2
%thetaU2_x_0=0; % not free to choose.
thetaU2_y_0=-pi/4;  
thetaU2_z_0=2*pi/3; % angle in the x-y plane
% thetaL2_x_0=0;
thetaL2_y_0=pi/5;
thetaL2_z_0=thetaU2_z_0;
%arm3
%thetaU3_x_0=0; % not free to choose.
thetaU3_y_0=-pi/4;  
thetaU3_z_0=0; % angle in the x-y plane
% thetaL3_x_0=0;
thetaL3_y_0=pi/5;
thetaL3_z_0=thetaU3_z_0; 

thetaA_0=[thetaU1_y_0 ; thetaU2_y_0 ;thetaU3_y_0]; % this is useful for Simulink

% -L_b/2-L_u*sin(thetaA1)+L_l*sin(thetaL_y_0)+L_p/2=0
L(2)=(L(4)/2-L(3)/2+L(1)*abs(sin(thetaU1_y_0)))/sin(thetaL1_y_0)
if L(2)<=0
    error('Invalid initial configuration. L(2)<0')
end
L_u=L(1); L_l=L(2); L_e=L(3); L_b=L(4);    

% Define rotation matrices
%arm1
R_A1=rot3D_Rodrigues([0;0;1],thetaU1_z_0); % constant
R_L1_0=R_A1*rot3D_Rodrigues([0;1;0],thetaL1_y_0); % constant
R_U1_0=R_A1*rot3D_Rodrigues([0;1;0],thetaA_0(1));
%arm2
R_A2=rot3D_Rodrigues([0;0;1],thetaU2_z_0); % constant
R_L2_0=R_A2*rot3D_Rodrigues([0;1;0],thetaL2_y_0); % constant
R_U2_0=R_A2*rot3D_Rodrigues([0;1;0],thetaA_0(2));
%arm3
R_A3=rot3D_Rodrigues([0;0;1],thetaU3_z_0); % constant
R_L3_0=R_A3*rot3D_Rodrigues([0;1;0],thetaL3_y_0); % constant
R_U3_0=R_A3*rot3D_Rodrigues([0;1;0],thetaA_0(3)); % symbolic

q0(1:3)=zeros(3,1) ; % intial M
% arm 1                                   
q0(4:6)= q0(1:3)+R_A1*[L(4)/2; 0 ;0]-R_U1_0*[0;0;L(1)/2];  % initial U1
q0(7:9)= [0;thetaU1_y_0;thetaU1_z_0];
q0(10:12)=q0(4:6) +R_U1_0*[0;0;-L(1)/2]-R_L1_0*[0;0;L(2)/2]; % initial L1
q0(13:15)=[0;thetaL1_y_0;thetaL1_z_0];
% arm 2
q0(16:18)= q0(1:3)+R_A2*[L(4)/2; 0 ;0]-R_U2_0*[0;0;L(1)/2];  % initial U2
q0(19:21)= [0;thetaU2_y_0;thetaU2_z_0];
q0(22:24)=q0(16:18) +R_U2_0*[0;0;-L(1)/2]-R_L2_0*[0;0;L(2)/2]; % initial L2
q0(25:27)=[0;thetaL2_y_0;thetaL2_z_0];
% arm 3
q0(28:30)= q0(1:3)+R_A3*[L(4)/2; 0 ;0]-R_U3_0*[0;0;L(1)/2];  % initial U3
q0(31:33)= [0;thetaU3_y_0;thetaU3_z_0];
q0(34:36)=q0(28:30) +R_U3_0*[0;0;-L(1)/2]-R_L3_0*[0;0;L(2)/2]; % initial L3
q0(37:39)=[0;thetaL3_y_0;thetaL3_z_0];
q0(40:42)=q0(34:36)+R_L3_0*[0;0;-L(2)/2]- R_A3*[L(3)/2;0 ;0]; % inital P
% Check P: the following should all be equal
% P1=q0(10:12)+R_L1_0*[0;0;-L(2)/2]- R_A1*[L(3)/2;0 ;0]
% P2=q0(22:24)+R_L2_0*[0;0;-L(2)/2]- R_A2*[L(3)/2;0 ;0]
% P3=q0(34:36)+R_L3_0*[0;0;-L(2)/2]- R_A3*[L(3)/2;0 ;0]

R_L0=zeros(3,3,3);
R_L0(:,:,1)=R_L1_0;
R_L0(:,:,2)=R_L2_0;
R_L0(:,:,3)=R_L3_0;
figure;
plot_Delta3D( q0,L,R_A1,R_A2,R_A3,R_L1_0,R_L2_0,R_L3_0,thetaA_0,0,0 )

qd0=q0(ind_d);

plot_boundaries(L,thetaU1_z_0,thetaU2_z_0,thetaU3_z_0);

%% Calculate and show an extended position
thetaA1=0.3;
P0=[(0.5*(L_b-L_e)-(L_u+L_l)*sin(thetaA1))*cos(thetaU1_z_0);
    (0.5*(L_b-L_e)-(L_u+L_l)*sin(thetaA1))*sin(thetaU1_z_0);
                  -(L_u+L_l)*cos(thetaA1)];
P0=0.999*P0;  
%P0=q0(40:42);
qa0=[thetaA1; thetaA_0(2);thetaA_0(3)]; 
[qa, ~, ~]=calculate_IK(L,qa0,P0,[0;0;0],[0;0;0],tolerance);

R_L0(:,:,1)=R_A1*rot3D_Rodrigues([0;1;0],thetaA1);
time_range=0;
[Qd,R1,R2,R3] = calculate_position_dep(q0,qa,R_L0,time_range,tolerance,L);

% Recover independent co-ordinates
Qi(1:3,:)=[0;thetaU1_y_0;thetaU1_z_0]+[-sin(thetaU1_z_0);cos(thetaU1_z_0);0]*(qa(1)-thetaA_0(1));
Qi(4:6,:)=[0;thetaU2_y_0;thetaU2_z_0]+[-sin(thetaU2_z_0);cos(thetaU2_z_0);0]*(qa(2)-thetaA_0(2));
Qi(7:9,:)=[0;thetaU3_y_0;thetaU3_z_0]+[-sin(thetaU3_z_0);cos(thetaU3_z_0);0]*(qa(3)-thetaA_0(3));

Q(ind_i,:)= Qi;
Q(ind_d,:)= Qd;

figure;
plot_Delta3D( Q,L,R_A1,R_A2,R_A3,R1,R2,R3,qa0,0,0 )
plot_boundaries(L,thetaU1_z_0,thetaU2_z_0,thetaU3_z_0);

