function Qacc_d=calculate_acceleration_dependent(Qacc_i,Qvel,R_L1_t,R_L2_t,R_L3_t,thetaA_t,L)
%Lior Sinai, 2019-02-23
% Caluculates accelerations for the dependent co-ordinates given known
% accelerations for the independent co-ordiantes and all velocities and
% positions
% INPUTS
% Qacc_i =  9xN independent co-ordinate accelerations
% Qvel   = 42xN velocites
% R_L1_t = 3x3xN rotation matrix for lower arm 1
% R_L2_t = 3x3xN rotation matrix for lower arm 2
% R_L3_t = 3x3xN rotation matrix for lower arm 3
% thetaA_t = 3xN actuator angles
% L=[L_upper L_lower L_endEffector L_base] ... lengths

% OUTPUTS
% Qacc_d = 33xN dependent co-ordinate accelerations

SIZE_QD=33;
SIZE_QI=9;
SIZE_Q=42;


N=size(Qacc_i,2);
Qacc_d=zeros(SIZE_QD,N);

% Extract values
L_u=L(1);
L_l=L(2);

fprintf('%d time steps. Progess of accelerations: 000.0%%\n',N)
for timeStep=1:N
    fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
    fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
    %Match symbols
    
    %match symbols to numeric values
    thetaA1=thetaA_t(1,timeStep);
    thetaA2=thetaA_t(2,timeStep);
    thetaA3=thetaA_t(3,timeStep);
    
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
    
    % calculate the Jacobians
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
   
   % calculate gamma
    qvel=Qvel(:,timeStep);

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

   gamma = Gamma_without_drivingNumeric(...
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
   
   % finally, calculate the dependent accelerations
   Qacc_d(:,timeStep)=-Jacobian_d\Jacobian_i*Qacc_i(:,timeStep)+...
                      +Jacobian_d\gamma;
   
                           
end