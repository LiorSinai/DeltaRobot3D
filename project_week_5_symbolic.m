%% Project Week 5: Equations of Motion

clear
% 18 January 2019
% Calculate all values symbolically

syms L_u L_l L_e L_b
L=[L_u L_l L_e L_b];  % L_upper, L_lower, L_effector, L_base
omega = sym('omega_',[3,1]);

% symbolic variables used
syms M_x U1_x U2_x U3_x L1_x L2_x L3_x P_x 
syms M_y U1_y U2_y U3_y L1_y L2_y L3_y P_y 
syms M_z U1_z U2_z U3_z L1_z L2_z L3_z P_z 
syms thetaU1_x thetaU2_x thetaU3_x 
syms thetaU1_y thetaU2_y thetaU3_y % constant
syms thetaU1_z thetaU2_z thetaU3_z % constant
syms thetaA1 thetaA2 thetaA3 % actuator angles used to parameterise matrix, not necessarily the same as the previous inertial frame angles
syms thetaL1_x thetaL2_x thetaL3_x
syms thetaL1_y thetaL2_y thetaL3_y
syms thetaL1_z thetaL2_z thetaL3_z
% lower angles are free, therefore the rotation matrix actually requires 4
% parameters: 3 co-ords for axis of rotation, and 1 for the actual rotation
%          or 3 angles and 1 for axis of rotation (implicit with Euler
%          angles, but may cause singularities -> crash program)
% So treat it as a black box by making it all symbols:
R_L1=sym('rL1_',[3 3]); 
R_L2=sym('rL2_',[3 3]);
R_L3=sym('rL3_',[3 3]);
syms t_sym

% Order is important. See equation 5.33 on pg 60
Q_sym=[M_x;  M_y;  M_z; 
    % arm 1
   U1_x; U1_y; U1_z; thetaU1_x; thetaU1_y; thetaU1_z;
   L1_x; L1_y; L1_z; thetaL1_x; thetaL1_y; thetaL1_z;
    % arm 2
   U2_x; U2_y; U2_z; thetaU2_x; thetaU2_y; thetaU2_z;
   L2_x; L2_y; L2_z; thetaL2_x; thetaL2_y; thetaL2_z;
    % arm 3
   U3_x; U3_y; U3_z; thetaU3_x; thetaU3_y; thetaU3_z;
   L3_x; L3_y; L3_z; thetaL3_x; thetaL3_y; thetaL3_z;
   P_x;  P_y; P_z
  ];
Qvel_sym = sym('qvel',[length(Q_sym),1]);

% Set initial rotation matrices
%%%             z
%%%             |
%%% Start at 1__M__2_ _ _x          y
%%%          /     \U           3\ |
%%%          \     /L             \|_ _2 x
%%%           \_P_/               /
%%%                             1/
%arm1
%thetaU1_x_0=0;
thetaU1_y_0=-pi/4;  
thetaU1_z_0=-2*pi/3; % angle in the x-y plane. Do not change
thetaL1_x_0=0;
thetaL1_y_0=pi/4;
thetaL1_z_0=thetaU1_z_0;
%arm2
%thetaU2_x_0=0;
thetaU2_y_0=-pi/4;  
thetaU2_z_0=2*pi/3; % angle in the x-y plane. Do not change
thetaL2_x_0=0;
thetaL2_y_0=pi/4;
thetaL2_z_0=thetaU2_z_0;
%arm3
%thetaU3_x_0=0;
thetaU3_y_0=-pi/4;  
thetaU3_z_0=0; % angle in the x-y plane. Do not change
thetaL3_x_0=0;
thetaL3_y_0=pi/4;
thetaL3_z_0=thetaU3_z_0; 

% Define rotation matrices
%arm1
R_A1=rot3D_Rodrigues([0;0;1],thetaU1_z_0); % constant
R_L1_0=R_A1*rot3D_Rodrigues([0;1;0],thetaL1_y_0); % constant
R_U1=R_A1*rot3D_Rodrigues([0;1;0],thetaA1); % symbolic
R_U1_0=eval(subs(R_U1,thetaA1,thetaU1_y_0));% constant
%arm2
R_A2=rot3D_Rodrigues([0;0;1],thetaU2_z_0); % constant
R_L2_0=R_A2*rot3D_Rodrigues([0;1;0],thetaL2_y_0); % constant
R_U2=R_A2*rot3D_Rodrigues([0;1;0],thetaA2); % symbolic
R_U2_0=eval(subs(R_U2,thetaA2,thetaU2_y_0));% constant
%arm3
R_A3=rot3D_Rodrigues([0;0;1],thetaU3_z_0); % constant
R_L3_0=R_A3*rot3D_Rodrigues([0;1;0],thetaL3_y_0); % constant
R_U3=R_A3*rot3D_Rodrigues([0;1;0],thetaA3); % symbolic
R_U3_0=eval(subs(R_U3,thetaA3,thetaU3_y_0));% constant

%% Solve Kinematics

Phi= [...
    ... # Arm 1:  3x3 = 9 constraints
    [M_x;M_y;M_z] + R_A1*[L(4)/2; 0 ;0] - [U1_x; U1_y; U1_z] - R_U1*[0;0;L(1)/2];
    [U1_x; U1_y; U1_z] + R_U1*[0;0;-L(1)/2] - [L1_x; L1_y; L1_z] - R_L1*[0;0;L(2)/2];
    [L1_x; L1_y; L1_z] + R_L1*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A1*[L(3)/2;0 ;0];
    ... # upper angles: 3 constraints
    % Important: R_U1(1,2) and R_U2(2,2) must not be functions of thetaU1
    % Therefore phi_0=thetaU1_x_0=0
    %     omega_U = R_U1*[0;omega_U;0]
    % so: theta_U = integral(R_U1*[0;omega_U;0])
    thetaU1_x - R_U1(1,2)*omega(1)*t_sym - 0; %thetaU1_x_0=0
    thetaU1_y - R_U1(2,2)*omega(1)*t_sym - thetaU1_y_0;% driving constraint
    thetaU1_z - 0 - thetaU1_z_0; % note: R(3,2)=0 if phi_0=thetaU1_x_0=0
    ... # Arm 2:  3x3 = 9 constraints
    [M_x;M_y;M_z] + R_A2*[L(4)/2; 0 ;0] - [U2_x; U2_y; U2_z] - R_U2*[0;0;L(1)/2];
    [U2_x; U2_y; U2_z] + R_U2*[0;0;-L(1)/2] - [L2_x; L2_y; L2_z] - R_L2*[0;0;L(2)/2];
    [L2_x; L2_y; L2_z] + R_L2*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A2*[L(3)/2;0 ;0];
    ... # upper angles: 3 constraints
    thetaU2_x - R_U2(1,2)*omega(2)*t_sym - 0;  
    thetaU2_y - R_U2(2,2)*omega(2)*t_sym - thetaU2_y_0;% driving constraint
    thetaU2_z - 0 - thetaU2_z_0; % note: this variable is not used    
    ... # Arm 3:  3x3 = 9 constraints
    [M_x;M_y;M_z] + R_A3*[L(4)/2; 0 ;0] - [U3_x; U3_y; U3_z] - R_U3*[0;0;L(1)/2];
    [U3_x; U3_y; U3_z] + R_U3*[0;0;-L(1)/2] - [L3_x; L3_y; L3_z] - R_L3*[0;0;L(2)/2];
    [L3_x; L3_y; L3_z] + R_L3*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A3*[L(3)/2;0 ;0];
    ... # upper angles: 3 constraints
    thetaU3_x - R_U3(1,2)*omega(3)*t_sym - 0  ;
    thetaU3_y - R_U3(2,2)*omega(3)*t_sym - thetaU3_y_0;% driving constraint
    thetaU3_z - 0 - thetaU3_z_0; % note: this variable is not used    
        ... # Fixed variables:  3 constraints
    [M_x;M_y;M_z];
     ... # Fixed lower link spin in own axis
%      thetaL1_z-thetaL1_z_0;
%      thetaL2_z-thetaL2_z_0;
%      thetaL3_z-thetaL3_z_0 % correct but doesn't do anything -> still get
%      %rotations (don't subs these angles in anywhere)
    [0 1 0]*R_L1.'*R_L1_0*[1;0;0];  % L1_y has no component in L1_x_0
    %[1 0 0]*R_L1.'*R_L1_0*[0;1;0]; % L1_x has no component in L1_y_0 % equally
    %valid but can't have both ??
    [0 1 0]*R_L2.'*R_L2_0*[1;0;0]; % L2_y has no component in L2_x_0
    [0 1 0]*R_L3.'*R_L3_0*[1;0;0]; % L3_y has no component in L3_x_0
];

% Calculate the Jacobian using Matlab toolbox. But will need to
% fix the parts related to R and omega with equation 5.33 on pg 60
Jacobian=jacobian(Phi,Q_sym);
% Jacobian*Qdot
% colum  1:3 
% Qdot =[M_v;
%         4:6  7:9     10:12    13:15 
%        U1_v U1_omega L1_v L1_omega 
%         16:18 19:             :27
%        U2_v U2_omega L2_v L2_omega
%         28:                  39
%        U3_v U3_omega L3_v L3_omega
%       P_v ]
Jacobian(1:3,7:9)    =-R_U1*tilde([0;0;+L(1)/2]).'*R_U1.'; % check: =omegaTilde+omegaTilde.'=0
Jacobian(4:6,7:9)    =+R_U1*tilde([0;0;-L(1)/2]).'*R_U1.';
Jacobian(4:6,13:15)  =-R_L1*tilde([0;0;+L(2)/2]).'*R_L1.';
Jacobian(7:9,13:15)  =+R_L1*tilde([0;0;-L(2)/2]).'*R_L1.';

Jacobian(13:15,19:21) =-R_U2*tilde([0;0;+L(1)/2]).'*R_U2.';
Jacobian(16:18,19:21) =+R_U2*tilde([0;0;-L(1)/2]).'*R_U2.';
Jacobian(16:18,25:27) =-R_L2*tilde([0;0;+L(2)/2]).'*R_L2.';
Jacobian(19:21,25:27) =+R_L2*tilde([0;0;-L(2)/2]).'*R_L2.';

Jacobian(25:27,31:33) =-R_U3*tilde([0;0;+L(1)/2]).'*R_U3.';
Jacobian(28:30,31:33) =+R_U3*tilde([0;0;-L(1)/2]).'*R_U3.';
Jacobian(28:30,37:39) =-R_L3*tilde([0;0;+L(2)/2]).'*R_L3.';
Jacobian(31:33,37:39) =+R_L3*tilde([0;0;-L(2)/2]).'*R_L3.';

%lower link Jacobian
Jacobian(40,13:15) = [0 1 0]*R_L1.'*tilde(R_L1_0*[1;0;0]);
Jacobian(41,25:27) = [0 1 0]*R_L2.'*tilde(R_L2_0*[1;0;0]);
Jacobian(42,37:39) = [0 1 0]*R_L3.'*tilde(R_L3_0*[1;0;0]);

%% Analyse singularities
Qp = [P_x,P_y,P_z];

% independent co-ordinates  
ind_i=[7:9,19:21,31:33]';
Q_i = Q_sym(ind_i);
Qa = [thetaA1, thetaA2, thetaA3];
% ActuatorPhi=[Qi(1:3)-[0;thetaU1_y_0;thetaU1_z_0]-R_U1*([0;thetaA1;0]-[0;thetaU1_y_0;0]);
%              Qi(4:6)-[0;thetaU2_y_0;thetaU2_z_0]-R_U2*([0;thetaA2;0]-[0;thetaU2_y_0;0]);
%              Qi(7:9)-[0;thetaU3_y_0;thetaU3_z_0]-R_U3*([0;thetaA3;0]-[0;thetaU3_y_0;0]);
%              ];
ActuatorPhi=[Q_i(1)-R_U1(1,2)*(thetaA1-thetaU1_y_0);
             Q_i(4)-R_U2(1,2)*(thetaA2-thetaU2_y_0);
             Q_i(8)-thetaU3_y_0-R_U3(2,2)*(thetaA3-thetaU3_y_0)];
% ActuatorPhi=[Q_i(2)-thetaU1_y_0-R_U1(2,2)*(thetaA1-thetaU1_y_0);
%              Q_i(2)-thetaU1_y_0-R_U2(2,2)*(thetaA2-thetaU2_y_0);
%              Q_i(8)-thetaU3_y_0-R_U3(2,2)*(thetaA3-thetaU3_y_0)];
             
% dependent co-ordinates
% Order is important. See equation 5.33 on pg 60
ind_d=(1:42)';
ind_d(ind_i)=[]; % remove independent co-ordinates
Q_d=Q_sym(ind_d);

Phi_without_driving=[Phi(1:9,:) ; Phi(13:21,:) ; Phi(25:33,:); Phi(37:42,:)];
Jacobian_without_driving=[Jacobian(1:9,:) ; Jacobian(13:21,:) ; Jacobian(25:33,:); Jacobian(37:42,:)];

Jacobian_d=Jacobian_without_driving(:,ind_d);
Jacobian_i=Jacobian_without_driving(:,ind_i);

Jacobian_a=-jacobian(ActuatorPhi,Qa)\jacobian(ActuatorPhi,Q_i);
% Jacobian_a =[ 1/(+3^(1/2)/2), 1/(-0.5), 0 0, 0, 0, 0, 0, 0;
%                0, 0, 0, 1/(-3^(1/2)/2),1/(-0.5), 0, 0, 0, 0;
%                0, 0, 0,                    0, 0, 0, 0, 1, 0];
rank(Jacobian_a);

%really slow, so load from workspace
% J=-Jacobian_d\Jacobian_i;
% Jp=J(31:33,:);
% J_FK=Jp*pinv(Jacobian_a); % jacobian forward kinematics. Need to take
% pseudo inverse since inverse does not exist. Also two solutions exist,
% since can calculate ActuatorPhi in two ways
% detJFK=det(J_FK);
% detJFK=simplify(detJFK);
% save('detJFK','detJFK');
load('detJFK');

% solution=solve(det(Jp),Qi); % crashes Matlab
% Solutions from 2D case:
% thetaU1=thetaL1; 
% So the following should be zero:
subs(detJFK,R_L3,R_U3)

% Off plane angles: thetaA1 and thetaA2 not in the x-y-z plane
subs(detJFK,R_L1,R_U1)
subs(detJFK,R_L2,R_U2)

%% Symoblic Gamma
Gamma_1 = sym(zeros(size(Jacobian)));

Gamma_1(1:3,7:9)    =-tilde(Qvel_sym(7:9))*R_U1*tilde([0;0;+L(1)/2]).'*R_U1.';
Gamma_1(4:6,7:9)    =+tilde(Qvel_sym(7:9))*R_U1*tilde([0;0;-L(1)/2]).'*R_U1.';
Gamma_1(4:6,13:15)  =-tilde(Qvel_sym(13:15))*R_L1*tilde([0;0;+L(2)/2]).'*R_L1.';
Gamma_1(7:9,13:15)  =+tilde(Qvel_sym(13:15))*R_L1*tilde([0;0;-L(2)/2]).'*R_L1.';

Gamma_1(13:15,19:21) =-tilde(Qvel_sym(19:21))*R_U2*tilde([0;0;+L(1)/2]).'*R_U2.';
Gamma_1(16:18,19:21) =+tilde(Qvel_sym(19:21))*R_U2*tilde([0;0;-L(1)/2]).'*R_U2.';
Gamma_1(16:18,25:27) =-tilde(Qvel_sym(25:27))*R_L2*tilde([0;0;+L(2)/2]).'*R_L2.';
Gamma_1(19:21,25:27) =+tilde(Qvel_sym(25:27))*R_L2*tilde([0;0;-L(2)/2]).'*R_L2.';

Gamma_1(25:27,31:33) =-tilde(Qvel_sym(31:33))*R_U3*tilde([0;0;+L(1)/2]).'*R_U3.';
Gamma_1(28:30,31:33) =+tilde(Qvel_sym(31:33))*R_U3*tilde([0;0;-L(1)/2]).'*R_U3.';
Gamma_1(28:30,37:39) =-tilde(Qvel_sym(37:39))*R_L3*tilde([0;0;+L(2)/2]).'*R_L3.';
Gamma_1(31:33,37:39) =+tilde(Qvel_sym(37:39))*R_L3*tilde([0;0;-L(2)/2]).'*R_L3.';

%lower link Jacobian
Gamma_1(40,13:15) = [0 1 0]*(tilde(Qvel_sym(13:15))*R_L1).'*tilde(R_L1_0*[1;0;0]);
Gamma_1(41,25:27) = [0 1 0]*(tilde(Qvel_sym(25:27))*R_L2).'*tilde(R_L2_0*[1;0;0]);
Gamma_1(42,37:39) = [0 1 0]*(tilde(Qvel_sym(37:39))*R_L3).'*tilde(R_L3_0*[1;0;0]);

Gamma = -Gamma_1*Qvel_sym;

Gamma_without_driving=[Gamma(1:9,:) ; Gamma(13:21,:) ; Gamma(25:33,:); Gamma(37:42,:)];
%Position_without_driving=[Position(1:9,:) ; Position(13:21,:) ; Position(25:33,:); Position(37:42,:)];

%% Equations of motion
% Paramters
syms m_u m_l  m_p %Upper link mass, lower link mass, end effector mass
syms I_upper_xx I_upper_yy I_upper_zz
syms I_lower_xx I_lower_yy I_lower_zz
syms g % Gravitational acceleration

% Calculate moment inertias in the inertia frame
I_U1=R_U1*diag([I_upper_xx,I_upper_yy I_upper_zz])*R_U1.';
I_U2=R_U2*diag([I_upper_xx,I_upper_yy I_upper_zz])*R_U2.';
I_U3=R_U3*diag([I_upper_xx,I_upper_yy I_upper_zz])*R_U3.';

I_L1=R_L1*diag([I_lower_xx,I_lower_yy I_lower_zz])*R_L1.';
I_L2=R_L2*diag([I_lower_xx,I_lower_yy I_lower_zz])*R_L2.';
I_L3=R_L3*diag([I_lower_xx,I_lower_yy I_lower_zz])*R_L3.';

% Q=[M_x;  M_y;  M_z; 
%     % arm 1
%    U1_x; U1_y; U1_z; thetaU1_x; thetaU1_y; thetaU1_z;
%    L1_x; L1_y; L1_z; thetaL1_x; thetaL1_y; thetaL1_z;
%     % arm 2
%    U2_x; U2_y; U2_z; thetaU2_x; thetaU2_y; thetaU2_z;
%    L2_x; L2_y; L2_z; thetaL2_x; thetaL2_y; thetaL2_z;
%     % arm 3
%    U3_x; U3_y; U3_z; thetaU3_x; thetaU3_y; thetaU3_z;
%    L3_x; L3_y; L3_z; thetaL3_x; thetaL3_y; thetaL3_z;
%    P_x;  P_y; P_z
%   ];
M_0= diag([0 0 0 ...
        m_u*[1 1 1] I_upper_xx I_upper_yy I_upper_zz ...
        m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
        m_u*[1 1 1] I_upper_xx I_upper_yy I_upper_zz ...
        m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
        m_u*[1 1 1] I_upper_xx I_upper_yy I_upper_zz ...
        m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
        m_p*[1 1 1] ]); 
% correct moment of inertias
M_0(7:9,7:9)=I_U1;
M_0(13:15,13:15)=I_L1;    
M_0(19:21,19:21)=I_U2;  
M_0(25:27,25:27)=I_L2;  
M_0(31:33,31:33)=I_U3;
M_0(37:39,37:39)=I_L3;

% Partition into M*Qddot=[Mdd Mdi; Mid Mii]*[Qd_ddot Qi_ddot]
%Qi = [thetaU1_y;thetaU2_y];
Mii = [I_U1 zeros(3,6);
      zeros(3,3) I_U2  zeros(3,3);
      zeros(3,6) I_U3];
% Qd=[M_x;  M_y;  M_z; 
%     % arm 1
%    U1_x; U1_y; U1_z; thetaU1_x;           thetaU1_z;
%    L1_x; L1_y; L1_z; thetaL1_x; thetaL1_y; thetaL1_z;
%     % arm 2
%    U2_x; U2_y; U2_z; thetaU2_x;            thetaU2_z;
%    L2_x; L2_y; L2_z; thetaL2_x; thetaL2_y; thetaL2_z;
%     % arm 3
%    U3_x; U3_y; U3_z; thetaU3_x;            thetaU3_z;
%    L3_x; L3_y; L3_z; thetaL3_x; thetaL3_y; thetaL3_z;
%    P_x;  P_y; P_z
%   ];
Mdd = diag([0 0 0 ...
           m_u*[1 1 1]              ...
           m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
           m_u*[1 1 1]              ...
           m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
           m_u*[1 1 1]              ...
           m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
           m_p*[1 1 1] ]); 
% Correct Mdd
Mdd(10:12,10:12)=I_L1;
Mdd(19:21,19:21)=I_L2;
Mdd(28:30,28:30)=I_L3;

Mdi = sym(zeros(33,9));
Mid = sym(zeros(9,33));


Q_Ad=sym(zeros(33,1));
% Gravity
Q_Ad(6)=-m_u*g;
Q_Ad(9)=-m_l*g;
Q_Ad(15)=-m_u*g;
Q_Ad(18)=-m_l*g;
Q_Ad(24)=-m_u*g;
Q_Ad(27)=-m_l*g;
Q_Ad(33)=-m_p*g;


% Calculate gyroscopic forces in the inertia frame
C_U1=tilde(Qvel_sym(7:9))*  I_U1*Qvel_sym(7:9);
C_U2=tilde(Qvel_sym(19:21))*I_U2*Qvel_sym(19:21);
C_U3=tilde(Qvel_sym(31:33))*I_U3*Qvel_sym(31:33);

C_L1=tilde(Qvel_sym(13:15))*I_L1*Qvel_sym(13:15);
C_L2=tilde(Qvel_sym(25:27))*I_L2*Qvel_sym(25:27);
C_L3=tilde(Qvel_sym(37:39))*I_L3*Qvel_sym(37:39);

% dependent C, lower angles
Q_Ad(10:12)=-C_L1;  
Q_Ad(19:21)=-C_L2;
Q_Ad(28:30)=-C_L3;

% independent C
Q_Ai=sym(zeros(9,1));
Q_Ai(1:3)=-C_U1;
Q_Ai(4:6)=-C_U2;
Q_Ai(7:9)=-C_U3;

% % really slow, so load from workspace
% J = Jacobian_d\Jacobian_i;
% J=simplify(J);
% save('J','J')
% load('J')

% P2 = Jacobian_d\Gamma_without_driving;
% 
% Mhat = Mii - Mid*J - J.'*(Mdi - Mdd*J); 
%Qhat = QA_i-Mid*P2- P.'*[QAd - Mdd*P2];

%% Inverse Kinematics
Phi_IK=[
    ((0.5*(L_b-L_e)-L_u*sin(thetaA1))*cos(thetaU1_z_0)-P_x)^2+...
    ((0.5*(L_b-L_e)-L_u*sin(thetaA1))*sin(thetaU1_z_0)-P_y)^2+...
    (-L_u*cos(thetaA1)-P_z)^2-L_l^2;
    ((0.5*(L_b-L_e)-L_u*sin(thetaA2))*cos(thetaU2_z_0)-P_x)^2+...
    ((0.5*(L_b-L_e)-L_u*sin(thetaA2))*sin(thetaU2_z_0)-P_y)^2+...
    (-L_u*cos(thetaA2)-P_z)^2-L_l^2;
    ((0.5*(L_b-L_e)-L_u*sin(thetaA3))*cos(thetaU3_z_0)-P_x)^2+...
    ((0.5*(L_b-L_e)-L_u*sin(thetaA3))*sin(thetaU3_z_0)-P_y)^2+...
    (-L_u*cos(thetaA3)-P_z)^2-L_l^2
    ];
Jacobian_P=jacobian(Phi_IK,[P_x P_y P_z]);
Jacobian_qa=jacobian(Phi_IK,[thetaA1 thetaA2 thetaA3]);

syms Pvel_x Pvel_y Pvel_z omegaA1 omegaA2 omegaA3
syms Pacc_x Pacc_y Pacc_z
Pvel=[Pvel_x Pvel_y Pvel_z].';
Pacc=[Pacc_x Pacc_y Pacc_z].';
omegaA=[omegaA1 omegaA2 omegaA3].';
% The following terms are equal to Jacobian*Qvel, not just the Jacobian
Jacobian_PP=jacobian(Jacobian_P*Pvel,[P_x;P_y;P_z])*Pvel;
Jacobian_qaqa=jacobian(Jacobian_qa*omegaA,[thetaA1;thetaA2;thetaA3])*omegaA;
Jcross=jacobian(Jacobian_qa*omegaA,[P_x;P_y;P_z])*Pvel;
%J_cross=jacobian(Jacobian_P*Pvel,[thetaA1;thetaA2;thetaA3])*omegaA
Jcross2=sym(zeros(3,3));
Jcross2(:,1)=jacobian(Jacobian_qa(:,1),[P_x;P_y;P_z])*Pvel;
Jcross2(:,2)=jacobian(Jacobian_qa(:,2),[P_x;P_y;P_z])*Pvel;
Jcross2(:,3)=jacobian(Jacobian_qa(:,3),[P_x;P_y;P_z])*Pvel;
Jcross2=Jcross2*omegaA;
% Jcross2(:,1)=jacobian(Jacobian_P(:,1),[thetaA1;thetaA2;thetaA3])*omegaA;
% Jcross2(:,2)=jacobian(Jacobian_P(:,2),[thetaA1;thetaA2;thetaA3])*omegaA;
% Jcross2(:,3)=jacobian(Jacobian_P(:,3),[thetaA1;thetaA2;thetaA3])*omegaA;
% Jcross2=Jcross2*Pvel;
simplify(Jcross2-Jcross)
Gamma_IK=-(Jacobian_qaqa+Jacobian_PP+2*Jcross+Jacobian_P*Pacc);
%% Convert the functions to a numeric format
%Mfile=matlabFunction(Mhat,'File','MhatNumeric');
%Qfile=matlabFunction(Qhat,'File','QhatNumeric');
%Jacobian_file=matlabFunction(Jacobian,'File','Jacobian_Numeric');
%Jacobian_dfile=matlabFunction(Jacobian_d,'File','Jacobian_dNumeric');
%Jacobian_ifile=matlabFunction(Jacobian_i,'File','Jacobian_iNumeric');
%Jacobian_Pfile=matlabFunction(Jacobian_P,'File','Jacobian_PNumeric');
%Jacobian_qafile=matlabFunction(Jacobian_qa,'File','Jacobian_qaNumeric');
%GammaFile= matlabFunction(Gamma,'File','GammaNumeric');
%GammaFile_without_driving= matlabFunction(Gamma_without_driving,'File','Gamma_without_drivingNumeric');
%Gamma_IKFile= matlabFunction(Gamma_IK,'File','Gamma_IKNumeric');
