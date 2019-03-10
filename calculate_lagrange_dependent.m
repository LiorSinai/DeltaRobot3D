function [lambda,Q_Ad_out]=calculate_lagrange_dependent(...
    Qacc_i,Qacc_d,Qvel,...
    R_A1,R_A2,R_A3,...
    R_L1_t,R_L2_t,R_L3_t,thetaA_t,...
    L,I_C,m_u,m_l,m_p,g)
% Caculates Langrange multipliers for every timestep in Q for the dependent
% co-ordinate. 

% INPUTS
% Qacc_i =  9xN independet co-ordinate accelerations
% Qacc_d = 33xN dependet co-ordinate accelerations
%   Qvel = 42xN co-ordinate velocities
%   R_A1 = 3x3 rotation matrix for base arm 1 about the z-axis
%   R_A2 = 3x3 rotation matrix for base arm 2 about the z-axis
%   R_A3 = 3x3 rotation matrix for base arm 3 about the z-axis
% R_L1_t = 3x3xN rotation matrix for lower arm 1
% R_L2_t = 3x3xN rotation matrix for lower arm 2
% R_L3_t = 3x3xN rotation matrix for lower arm 3
% thetaA_t = 3xN actuator angle values
% L=[L_upper L_lower L_endEffector L_base] ... lengths [m]
%    I_C = 6x1 list of moment of inertias for the upper and lower arm
%    m_u = mass of the upper link [kg]
%    m_l = mass of the lower link [kg]
%    m_p = mass of the end effector (point P) [kg]
%     g = gravitational acceleration

% OUTPUTS
%  lambda = 33xN Lagrange multipliers at each time step
%    Q_AC = Q_A-Q_C = 33xN applied forces and coriolis forces matrix

SIZE_QD=33;

N=size(thetaA_t,2);

lambda=zeros(SIZE_QD,N);
Q_Ad_out=lambda;

% Extract values
L_l=L(2);

I_upper_xx=I_C(1);
I_upper_yy=I_C(2);
I_upper_zz=I_C(3);
I_lower_xx=I_C(4);
I_lower_yy=I_C(5);
I_lower_zz=I_C(6);

fprintf('%d time steps. Progess of lambdas: 000.0%%\n',N)
for timeStep = 1:N
    %fprintf('% d ... %.1f%% complete of lambdas\n',timeStep,100*timeStep/N)
    fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
    fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
    
    %match symbols to numeric values
    thetaA1=thetaA_t(1,timeStep);
    thetaA2=thetaA_t(2,timeStep);
    thetaA3=thetaA_t(3,timeStep);
    
   % upper rotation matrices
    R_U1=R_A1*rot3D_Rodrigues([0;1;0],thetaA1); 
    R_U2=R_A2*rot3D_Rodrigues([0;1;0],thetaA2); 
    R_U3=R_A3*rot3D_Rodrigues([0;1;0],thetaA3); 

    % R_L1
    R_L1=R_L1_t(:,:,timeStep);
    rL1_1_1 = R_L1(1,1); rL1_1_2 = R_L1(1,2); rL1_1_3 = R_L1(1,3);
    rL1_2_1 = R_L1(2,1); rL1_2_2 = R_L1(2,2); rL1_2_3 = R_L1(2,3);
    rL1_3_1 = R_L1(3,1); rL1_3_2 = R_L1(3,2); rL1_3_3 = R_L1(3,3);

    % R_L2
    R_L2=R_L2_t(:,:,timeStep);
    rL2_1_1 = R_L2(1,1); rL2_1_2 = R_L2(1,2); rL2_1_3 = R_L2(1,3);
    rL2_2_1 = R_L2(2,1); rL2_2_2 = R_L2(2,2); rL2_2_3 = R_L2(2,3);
    rL2_3_1 = R_L2(3,1); rL2_3_2 = R_L2(3,2); rL2_3_3 = R_L2(3,3);

    % R_L3
    R_L3=R_L3_t(:,:,timeStep);
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
    
    % calculate the dependent Jacobian
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
    %% Compute dependent forces Q_Ad
    Q_Ad=zeros(33,1); % reset Q_Ad every iteration
    % Gravity
    Q_Ad(6)=-m_u*g;
    Q_Ad(9)=-m_l*g;
    Q_Ad(15)=-m_u*g;
    Q_Ad(18)=-m_l*g;
    Q_Ad(24)=-m_u*g;
    Q_Ad(27)=-m_l*g;
    Q_Ad(33)=-m_p*g;

    % Calculate gyroscopic forces in the inertia frame
    C_L1=tilde(Qvel(13:15,timeStep))*I_L1*Qvel(13:15,timeStep);
    C_L2=tilde(Qvel(25:27,timeStep))*I_L2*Qvel(25:27,timeStep);
    C_L3=tilde(Qvel(37:39,timeStep))*I_L3*Qvel(37:39,timeStep);
    % dependent C, lower angles
    Q_Ad(10:12)=Q_Ad(10:12)-C_L1;  
    Q_Ad(19:21)=Q_Ad(19:21)-C_L2;
    Q_Ad(28:30)=Q_Ad(28:30)-C_L3;
                                  
    lambda(:,timeStep)=Jacobian_d.'\(Q_Ad-Mdd*Qacc_d(:,timeStep)-Mdi*Qacc_i(:,timeStep));
    
    Q_Ad_out(:,timeStep)=Q_Ad;

end