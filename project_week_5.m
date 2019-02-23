
clear
timestamp=datestr(now,'yyyymmdd_HHMMSS');
addpath('Numeric Functions')

% Control settings parameters
calculate = true;
plotBoolean = true;

% plot settings parameters
plot_positionBoolean = true;
plot_velocityBoolean = true;
plot_accelerationBoolean = true;

% input parameters
L = [0.5,1,0.3,0.8]; % L_upper, L_lower, L_effector, L_base

% Push down settings
% omega = 0.5*[1;1;1];
% delta_t = 0.1;
% delay_between_plots = delta_t;
% tolerance = 0.0001;
% end_time = 2;

% Twist settings
% omega = 0.5*[1;0.1;1];
% delta_t = 0.01;
% delay_between_plots = delta_t;
% tolerance = 0.00001;
% end_time = 3;

% Singularity settings
omega = 0.5*[1;0;0];
delta_t = 0.01;
delay_between_plots = delta_t;
tolerance = 0.0001;
end_time = 4;

% Full cycle settings
% L = [1,1,0.8,0.8]; % L_upper, L_lower, L_effector, L_base
% omega = 0.5*[1;0;0];
% delta_t = 0.01;
% delay_between_plots = delta_t;
% tolerance = 0.0001;
% end_time = 12.6;

% latent variables
if calculate == true
    time_range = 0:delta_t:end_time;
    Q = zeros(42,length(time_range));
end

%% symbolic variables used
syms M_x U1_x U2_x U3_x L1_x L2_x L3_x P_x 
syms M_y U1_y U2_y U3_y L1_y L2_y L3_y P_y 
syms M_z U1_z U2_z U3_z L1_z L2_z L3_z P_z 
syms thetaU1_x thetaU2_x thetaU3_x 
syms thetaU1_y thetaU2_y thetaU3_y % constant
syms thetaA1 thetaA2 thetaA3 % actuator angles used to parameterise matrix, not necessarily the same as the previous inertial frame angles
syms thetaU1_z thetaU2_z thetaU3_z % constant
syms thetaL1_x thetaL2_x thetaL3_x
syms thetaL1_y thetaL2_y thetaL3_y
syms thetaL1_z thetaL2_z thetaL3_z
% lower angles are free, therefore the rotation matrix actually requires 4
% parameters: 3 co-ords for axis of rotation, and 1 for the actual rotation
% So treat it as a black box by making it all symbols:
R_L1=sym('rL1_',[3 3]); 
R_L2=sym('rL2_',[3 3]);
R_L3=sym('rL3_',[3 3]);
syms t_sym

%syms M_x_dot U1_x_dot U2_x_dot L1_x_dot L2_x_dot P_x_dot M_y_dot U1_y_dot U2_y_dot L1_y_dot L2_y_dot P_y_dot thetaU1_dot thetaU2_dot thetaL1_dot thetaL2_dot

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

%% Set a good initial position
%%%             z
%%%             |
%%% Start at 1__M__2_ _ _x          y
%%%          /     \U           3\ |
%%%          \     /L             \|_ _2 x
%%%           \_P_/               /
%%%                             1/
q0=ones(size(Q_sym));
%arm1
%arm1
%thetaU1_x_0=0; not free to choose. This will cause problems if this is not
%zero, because then R(1,2) and R(2,2) will be a function of theta_U1_y
%See the driving constraints
thetaU1_y_0=-pi/4;  
thetaU1_z_0=-2*pi/3; % angle in the x-y plane
% thetaL1_x_0=0;
thetaL1_y_0=pi/4;
thetaL1_z_0=thetaU1_z_0;
%arm2
%thetaU2_x_0=0; % not free to choose.
thetaU2_y_0=-pi/4;  
thetaU2_z_0=2*pi/3; % angle in the x-y plane
% thetaL2_x_0=0;
thetaL2_y_0=pi/4;
thetaL2_z_0=thetaU2_z_0;
%arm3
%thetaU3_x_0=0; % not free to choose.
thetaU3_y_0=-pi/4;  
thetaU3_z_0=0; % angle in the x-y plane
% thetaL3_x_0=0;
thetaL3_y_0=pi/4;
thetaL3_z_0=thetaU3_z_0; 

% -L_b/2-L_u*sin(thetaA1)+L_l*sin(thetaL_y_0)+L_p/2=0
L(2)=(L(4)/2-L(3)/2+L(1)*abs(sin(thetaU1_y_0)))/sin(thetaL1_y_0)
if L(2)<=0
    error('Invalid initial configuration. L(2)<0')
end
    

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
% P1=q0(10:12)+R_L1_0*[0;0;-L(2)/2]- R_A1*[L(3)/2;0 ;0];
% P2=q0(22:24)+R_L2_0*[0;0;-L(2)/2]- R_A2*[L(3)/2;0 ;0];
% P3=q0(34:36)+R_L3_0*[0;0;-L(2)/2]- R_A3*[L(3)/2;0 ;0];

% independent and dependent co-ordinate indices
ind_i=[7:9,19:21,31:33]';
ind_d=(1:42)';
ind_d(ind_i)=[]; % remove independent co-ordinates

q0d = q0(ind_d);
q0i = q0(ind_i);

%% Rotation Matrices
R_sym=sym(zeros(3,3,6));
R_sym(:,:,1)=R_U1;
R_sym(:,:,2)=R_U2;
R_sym(:,:,3)=R_U3;
R_sym(:,:,4)=R_L1;
R_sym(:,:,5)=R_L2;
R_sym(:,:,6)=R_L3;

R_L0=zeros(3,3,3);
R_L0(:,:,1)=R_L1_0;
R_L0(:,:,2)=R_L2_0;
R_L0(:,:,3)=R_L3_0;

thetaA_sym=[thetaA1 ; thetaA2 ; thetaA3];
thetaA_0=[thetaU1_y_0 ; thetaU2_y_0 ;thetaU3_y_0]; % TODO: do these always have to be the same? If yes, then can remove this variable

plot_Delta3D( q0,L,R_A1,R_A2,R_A3,R_sym,thetaA_sym,R_L1_0,R_L2_0,R_L3_0,thetaA_0,0,0 )

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
    thetaU2_x - R_U2(1,2)*omega(2)*t_sym - 0; %thetaU2_x_0=0
    thetaU2_y - R_U2(2,2)*omega(2)*t_sym - thetaU2_y_0;% driving constraint
    thetaU2_z - 0 - thetaU2_z_0; % note: this variable is not used    
    ... # Arm 3:  3x3 = 9 constraints
    [M_x;M_y;M_z] + R_A3*[L(4)/2; 0 ;0] - [U3_x; U3_y; U3_z] - R_U3*[0;0;L(1)/2];
    [U3_x; U3_y; U3_z] + R_U3*[0;0;-L(1)/2] - [L3_x; L3_y; L3_z] - R_L3*[0;0;L(2)/2];
    [L3_x; L3_y; L3_z] + R_L3*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A3*[L(3)/2;0 ;0];
    ... # upper angles: 3 constraints
    thetaU3_x - R_U3(1,2)*omega(3)*t_sym - 0 ; %thetaU3_x_0=0
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
% Qvel =[M_v;
%         4:6  7:9     10:12    13:15 
%        U1_v U1_omega L1_v L1_omega 
%         16:18 19:             :27
%        U2_v U2_omega L2_v L2_omega
%         28:                  39
%        U3_v U3_omega L3_v L3_omega
%       P_v ]
Jacobian(1:3,7:9)    =-R_U1*tilde([0;0;+L(1)/2]).'*R_U1.';
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

Phi_t  = diff(Phi,t_sym); 
Phi_tt = diff(Phi_t,t_sym); 
Velocity = -Jacobian\Phi_t; %very slow line of code

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
% Acceleration = Jacobian\(-Phi_tt+Gamma);%very slow complex calculation. Not worth doing

%% Solve kinematic equations

if calculate == true
    fprintf('Starting calculations...\n')
    % Calculate for all time steps
%     tic
%     [Q,R1,R2,R3,thetaA_t] = calculate_position( Phi,Jacobian,Q_sym,R_sym,thetaA_sym,t_sym,q0,thetaA_0,R_L0,time_range,tolerance );
%     Qvel = calculate_velocity(Velocity,Q_sym,R_sym,thetaA_sym,t_sym,Q,R1,R2,R3,thetaA_t,time_range);    
%     QaccOld = calculate_acceleration(Jacobian,Phi_tt,Q_sym,R_sym,thetaA_sym,t_sym,Q,R1,R2,R3,Qvel,thetaA_t,time_range,L);
%     toc

    tic
    [Q,R1,R2,R3,thetaA_t] = calculate_position_fast(q0,thetaA_0,R_L0,time_range,tolerance,omega,L);
    Qvel = calculate_velocity_fast(Phi_t,t_sym,R1,R2,R3,thetaA_t,time_range,L);
    Qacc=calculate_acceleration_fast(Phi_tt,t_sym,R1,R2,R3,Qvel,thetaA_t,time_range,L);
    toc
    
   fprintf('Calculations complete\n')
end
%save(sprintf('TestRun_%s',timestamp)); % save all variables

%% Plot Matrix
if plotBoolean == true
   Nvalid=size(thetaA_t,2); % length of
   time_valid=time_range(1:Nvalid); % if a singularity was hit, the valid_time
                               % will be shorter than the time_range 
    if plot_velocityBoolean == true
        plot_velocity(Qvel,time_valid);
    end
    if plot_accelerationBoolean == true
        plot_accelerations2(Qacc,time_valid);
        plot_tangents(Qvel, Qacc,time_valid);
    end
     if plot_positionBoolean == true
        plot_Delta3D( Q,L,R_A1,R_A2,R_A3,R_sym,thetaA_sym,R1,R2,R3,thetaA_t,time_valid,delay_between_plots )
    end
end

%% Kinetic Parameters
m_u=0.35; % kg. Upper link mass
m_l=0.35; % kg. Lower link mass
I_upper_xx=(1/12)*m_u*L(1)^2;
I_upper_yy=(1/12)*m_u*L(1)^2;
I_upper_zz=I_upper_xx/10; 
I_lower_xx=(1/12)*m_l*L(2)^2;
I_lower_yy=(1/12)*m_l*L(2)^2;
I_lower_zz=I_lower_xx/10;
m_p   =1; % kg. end effector mass
g=9.81; % m/s^2. Gravitational acceleration

I_C=[I_upper_xx;I_upper_yy;I_upper_zz;...
     I_lower_xx;I_lower_yy;I_lower_zz];
M_C= diag([0 0 0 ...
        ... % arm 1
        m_u*[1 1 1] I_upper_xx I_upper_yy I_upper_zz ...
        m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
        ... % arm 2
        m_u*[1 1 1] I_upper_xx I_upper_yy I_upper_zz ...
        m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
        ... % arm 3
        m_u*[1 1 1] I_upper_xx I_upper_yy I_upper_zz ...
        m_l*[1 1 1] I_lower_xx I_lower_yy I_lower_zz ...
        ... % end effector
        m_p*[1 1 1] ]); 
    
Q_A=zeros(42,1); % applied forces
%Gravity
Q_A(6) =-m_u*g;
Q_A(12)=-m_l*g;
Q_A(18)=-m_u*g;
Q_A(24)=-m_l*g;
Q_A(30)=-m_u*g;
Q_A(36)=-m_l*g;
Q_A(42)=-m_p*g;

%% PID Controller
emax=0.1; % m = 0.1mm. maximum error
cycleTime = 0.35; %s for 1 kg
startDistance=0.025; %m
endDistance=0.305; %m
cycleDistance=endDistance-startDistance; %m
cycleAngle=atan2(cycleDistance,L(1)+L(2)); %=0.21*pi~pi/4

% reference constants
%hm=45*pi/180; % required angle
hm=thetaA_0(1)*0.5;
tm=cycleTime; % required time to reach hm
h0=thetaA_0(1);

% controller constansts
alpha=0.15;
beta=2;

% calculate cross over frequency
wc=(4*pi^2*beta*abs(hm)/(alpha*emax))^(1/3)/tm; % rad/s 
wc_rpm=wc * (60/(2*pi) ); % rpm

M =1 ; % there will be a multiplication by a full M matrix in the Simulink model

kc=M*wc^2*sqrt(alpha);
tz=1/(wc*sqrt(alpha));
tp=alpha*tz;
ti=beta*tz;

