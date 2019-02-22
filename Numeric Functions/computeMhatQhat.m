function [Mhat, Qhat,qd_vel,Gamma_without_driving,qvel] = ...
    computeMhatQhat(qi,qi_vel,R_L1,R_L2,R_L3,...
        m_u, m_p, m_l, L, g,...
        I_upper_xx, I_upper_yy, I_upper_zz,...
        I_lower_xx, I_lower_yy, I_lower_zz,...
        thetaA_0)
% Compute dynamics for the system 
% Assume the model is perfect, so this is a copy past of M_Q_AC_model

%Pre-initialise for Simulink
Mhat=zeros(9,9);
Qhat=zeros(9,1);
qd_vel=zeros(33,1);

% independent and dependent co-ordinate indices
ind_i=[7:9,19:21,31:33]';
% ind_d=(1:42)'; 
% ind_d(ind_i)=[]; % remove independent co-ordinates -> doesn't work in Simulink
ind_d = [1:6,10:18,22:30,34:42]';

%% Parameters
L_u=L(1);
L_l=L(2);
L_e=L(3);
L_b=L(4);

thetaA1 = qi(1)/(+3^(1/2)/2) + thetaA_0(1); 
thetaA2 = qi(4)/(-3^(1/2)/2) + thetaA_0(2); 
thetaA3 = qi(8); 
   
% upper rotation matrices
R_A1=rot3D_Rodrigues([0;0;1],-2*pi/3); % constant
R_U1=R_A1*rot3D_Rodrigues([0;1;0],thetaA1); 

R_A2=rot3D_Rodrigues([0;0;1],2*pi/3); % constant
R_U2=R_A2*rot3D_Rodrigues([0;1;0],thetaA2); 

R_A3=rot3D_Rodrigues([0;0;1],0); % constant
R_U3=R_A3*rot3D_Rodrigues([0;1;0],thetaA3); 

% R_L1
rL1_1_1 = R_L1(1,1); rL1_1_2 = R_L1(1,2); rL1_1_3 = R_L1(1,3);
rL1_2_1 = R_L1(2,1); rL1_2_2 = R_L1(2,2); rL1_2_3 = R_L1(2,3);
rL1_3_1 = R_L1(3,1); rL1_3_2 = R_L1(3,2); rL1_3_3 = R_L1(3,3);

% R_L2
rL2_1_1 = R_L2(1,1); rL2_1_2 = R_L2(1,2); rL2_1_3 = R_L2(1,3);
rL2_2_1 = R_L2(2,1); rL2_2_2 = R_L2(2,2); rL2_2_3 = R_L2(2,3);
rL2_3_1 = R_L2(3,1); rL2_3_2 = R_L2(3,2); rL2_3_3 = R_L2(3,3);

% R_L3
rL3_1_1 = R_L3(1,1); rL3_1_2 = R_L3(1,2); rL3_1_3 = R_L3(1,3);
rL3_2_1 = R_L3(2,1); rL3_2_2 = R_L3(2,2); rL3_2_3 = R_L3(2,3);
rL3_3_1 = R_L3(3,1); rL3_3_2 = R_L3(3,2); rL3_3_3 = R_L3(3,3);

%% Calculate moment inertias in the inertia frame
I_U1=R_U1*diag([I_upper_xx,I_upper_yy I_upper_zz])*R_U1.';
I_U2=R_U2*diag([I_upper_xx,I_upper_yy I_upper_zz])*R_U2.';
I_U3=R_U3*diag([I_upper_xx,I_upper_yy I_upper_zz])*R_U3.';

I_L1=R_L1*diag([I_lower_xx,I_lower_yy I_lower_zz])*R_L1.';
I_L2=R_L2*diag([I_lower_xx,I_lower_yy I_lower_zz])*R_L2.';
I_L3=R_L3*diag([I_lower_xx,I_lower_yy I_lower_zz])*R_L3.';

%% Partition into M*Qddot=[Mdd Mdi; Mid Mii]*[Qi_ddot Qd_ddot]
Mii = [I_U1 zeros(3,6);
      zeros(3,3) I_U2  zeros(3,3);
      zeros(3,6) I_U3];
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
  
Mdi = zeros(33,9);
Mid = zeros(9,33);

Jacobian_d = Jacobian_dNumeric(L_l,...
                               rL1_1_1,rL1_1_2,...%R_L1
                               rL1_2_1,rL1_2_2,...
                               rL1_3_1,rL1_3_2,...
                               rL2_1_1,rL2_1_2,...%R_L2
                               rL2_2_1,rL2_2_2,...
                               rL2_3_1,rL2_3_2,...
                               rL3_1_1,rL3_1_2,...%R_L3
                               rL3_2_1,rL3_2_2,...
                               rL3_3_1,rL3_3_2);
Jacobian_i = Jacobian_iNumeric(L_u,thetaA1,thetaA2,thetaA3);

%% Compute dependent velocities
J = Jacobian_d\Jacobian_i;
qd_vel = -J*qi_vel;

%match symbolic variables to numeric variables
qvel=zeros(42,1);
qvel(ind_i)=qi_vel;
qvel(ind_d)=qd_vel;

qvel1=qvel(1);qvel2=qvel(2);qvel3=qvel(3);
qvel4=qvel(4);qvel5=qvel(5);qvel6=qvel(6);
qvel7=qvel(7);qvel8=qvel(8);qvel9=qvel(9);
qvel10=qvel(10);qvel11=qvel(11);qvel12=qvel(12);
qvel13=qvel(13);qvel14=qvel(14);qvel15=qvel(15);
qvel16=qvel(16);qvel17=qvel(17);qvel18=qvel(18);
qvel19=qvel(19);qvel20=qvel(20);qvel21=qvel(21);
qvel22=qvel(22);qvel23=qvel(23);qvel24=qvel(24);
qvel25=qvel(25);qvel26=qvel(26);qvel27=qvel(27);
qvel28=qvel(28);qvel29=qvel(29);qvel30=qvel(30);
qvel31=qvel(31);qvel32=qvel(32);qvel33=qvel(33);
qvel34=qvel(34);qvel35=qvel(35);qvel36=qvel(36);
qvel37=qvel(37);qvel38=qvel(38);qvel39=qvel(39);
qvel40=qvel(40);qvel41=qvel(41);qvel42=qvel(42);

Gamma_without_driving = Gamma_without_drivingNumeric(...
    L_l,L_u,...
    qvel7,qvel8,qvel9,...
    qvel13,qvel14,qvel15,...
    qvel19,qvel20,qvel21,...
    qvel25,qvel26,qvel27,...
    qvel31,qvel32,qvel33,...
    qvel37,qvel38,qvel39,...
    rL1_1_1,rL1_1_2,...%R_L1
    rL1_2_1,rL1_2_2,...
    rL1_3_1,rL1_3_2,...
    rL2_1_1,rL2_1_2,...%R_L2
    rL2_2_1,rL2_2_2,...
    rL2_3_1,rL2_3_2,...
    rL3_1_1,rL3_1_2,...%R_L3
    rL3_2_1,rL3_2_2,...
    rL3_3_1,rL3_3_2,...
    thetaA1,thetaA2,thetaA3);

%% Compute Q_Ai and Q_Ad
Q_Ad=zeros(33,1);
% Gravity
Q_Ad(6)=-m_u*g;
Q_Ad(9)=-m_l*g;
Q_Ad(15)=-m_u*g;
Q_Ad(18)=-m_l*g;
Q_Ad(24)=-m_u*g;
Q_Ad(27)=-m_l*g;
Q_Ad(33)=-m_p*g;

% Calculate gyroscopic forces in the inertia frame
C_U1=tilde(qvel(7:9))*  I_U1*qvel(7:9);
C_U2=tilde(qvel(19:21))*I_U2*qvel(19:21);
C_U3=tilde(qvel(31:33))*I_U3*qvel(31:33);

C_L1=tilde(qvel(13:15))*I_L1*qvel(13:15);
C_L2=tilde(qvel(25:27))*I_L2*qvel(25:27);
C_L3=tilde(qvel(37:39))*I_L3*qvel(37:39);
% dependent C, lower angles
Q_Ad(10:12)=Q_Ad(10:12)-C_L1;  
Q_Ad(19:21)=Q_Ad(19:21)-C_L2;
Q_Ad(28:30)=Q_Ad(28:30)-C_L3;

% independent C
Q_Ai=zeros(9,1);
Q_Ai(1:3)=Q_Ai(1:3)-C_U1;
Q_Ai(4:6)=Q_Ai(4:6)-C_U2;
Q_Ai(7:9)=Q_Ai(7:9)-C_U3;
%% Compute Ouputs
J2 = Jacobian_d\Gamma_without_driving;

Mhat = Mii - Mid*J - J.'*(Mdi - Mdd*J); 
Qhat = Q_Ai-Mid*J2 -(J.')*(Q_Ad -Mdd*J2);

%M_0=eye(9,9);

end

