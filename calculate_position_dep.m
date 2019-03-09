function [Qd,R1_out,R2_out,R3_out] = calculate_position_dep(q0,thetaA_t,R_L0,time_range,tolerance,L)
% calulation for the dependent position co-ordinates
% does not use the Symbloic toolbox (but that is used to make the Jacobian function)
% This is less general than the normal code, as it uses a pre-defined system matrix
% and Jacobian

% INPUTS
%  q0 = 42x1 initial position (guess). This is converted to q0d in the code
% thetaA_t = 3xN actuator angle values
%  R_L0 = (3x3)x3 set of 3x3 rotation matrices. Initial values (guess)
% time_range = 1xN time values. 
% tolerance = desired accuracy for the Newton-Raphson iterations
% L=[L_upper L_lower L_endEffector L_base] ... lengths [m]

% OUTPUTS
% Note: if a singularity is hit, N is changed to N=singularity time step
% Qd = 33xN co-ordinate values
% R1_out= 3x3xN rotation matrices set for lower arm 1
% R2_out= 3x3xN rotation matrices set for lower arm 2
% R3_out= 3x3xN rotation matrices set for lower arm 3

%% initialise variables
SIZE_QD=33;
Qd=zeros(SIZE_QD,length(time_range));
R1_out=zeros(3,3,length(time_range));
R2_out=zeros(3,3,length(time_range));
R3_out=zeros(3,3,length(time_range));
singularity=false;
maxIteration=20;


%% extract initial values
R_L1_0=R_L0(:,:,1);
R_L2_0=R_L0(:,:,2);
R_L3_0=R_L0(:,:,3);

thetaU1_y_0=q0(8);  thetaU1_z_0=q0(9);
thetaU2_y_0=q0(20); thetaU2_z_0=q0(21);
thetaU3_y_0=q0(32); thetaU3_z_0=q0(33);

% independent and dependent co-ordinate indices
ind_i=[7:9,19:21,31:33]';
ind_d=(1:42)';
ind_d(ind_i)=[]; % remove independent co-ordinates
qd0=q0(ind_d);

% recalculate R_A matrices
R_A1=rot3D_Rodrigues([0;0;1],thetaU1_z_0); % constant
R_A2=rot3D_Rodrigues([0;0;1],thetaU2_z_0); % constant
R_A3=rot3D_Rodrigues([0;0;1],thetaU3_z_0); % constant

%% Newton Raphson Calculations
N=length(time_range);
fprintf('%d time steps. Progess of positions: 000.0%%\n',N)
for k = 1:N
    if ~singularity 
        %fprintf('% d ... %.1f%% complete of positions\n',k,100*k/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',k/N*100); %write the new number
    end
    if k == 1
       qd=qd0;
       R_L1 = R_L1_0; 
       R_L2 = R_L2_0; 
       R_L3 = R_L3_0; 
    end
    % set the error>tolerance on the first loop
    S =ones(SIZE_QD)*100*tolerance;
    newton_iterations = 0;

    while (norm(S) > tolerance)&& (~singularity)
        newton_iterations = newton_iterations +1;
        if newton_iterations>maxIteration
            warning('iterions = %d. Probably hit a singularity\naborting calculations\n',newton_iterations)
            % shorten time_range
            time_range=time_range(1:k);
            Qd=Qd(:,1:k);
            R1_out=R1_out(:,:,1:k);
            R2_out=R2_out(:,:,1:k);
            R3_out=R3_out(:,:,1:k);
            singularity=true;
            break;
        end
        % Match symbols to numeric values
        M_x = qd(1);  M_y = qd(2);  M_z = qd(3);
        U1_x = qd(4); U1_y = qd(5); U1_z = qd(6);
        %thetaU1_x=qd(7);thetaU1_y=qd(8);thetaU1_z=qd(9);
        L1_x = qd(7); L1_y = qd(8); L1_z = qd(9);
        thetaL1_x = qd(10); thetaL1_y = qd(11); thetaL1_z = qd(12);
        U2_x = qd(13); U2_y = qd(14); U2_z = qd(15);
        %thetaU2_x=qd(19);thetaU2_y=qd(20);thetaU2_z=qd(21);
        L2_x = qd(16); L2_y = qd(17); L2_z = qd(18);
        thetaL2_x = qd(19); thetaL2_y = qd(20); thetaL2_z = qd(21);
        U3_x = qd(22); U3_y = qd(23); U3_z = qd(24);
        %thetaU3_x=qd(31);thetaU3_y=qd(32);thetaU3_z=qd(33);
        L3_x = qd(25); L3_y = qd(26); L3_z = qd(27);
        thetaL3_x = qd(28); thetaL3_y = qd(29); thetaL3_z = qd(30);
        P_x = qd(31);   P_y = qd(32);  P_z = qd(33);
        t_sym=time_range(k);

        thetaA1 = thetaA_t(1,k);
        thetaA2 = thetaA_t(2,k);
        thetaA3 = thetaA_t(3,k);

        R_U1 =R_A1*rot3D_Rodrigues([0;1;0],thetaA1);         
        R_U2 =R_A2*rot3D_Rodrigues([0;1;0],thetaA2); 
        R_U3 =R_A3*rot3D_Rodrigues([0;1;0],thetaA3); 
        
        % Calculate the position matrix
        S = [...
        ... # Arm 1:  3x3 = 9 constraints
        [M_x;M_y;M_z] + R_A1*[L(4)/2; 0 ;0] - [U1_x; U1_y; U1_z] - R_U1*[0;0;L(1)/2];
        [U1_x; U1_y; U1_z] + R_U1*[0;0;-L(1)/2] - [L1_x; L1_y; L1_z] - R_L1*[0;0;L(2)/2];
        [L1_x; L1_y; L1_z] + R_L1*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A1*[L(3)/2;0 ;0];
        ... # Arm 2:  3x3 = 9 constraints
        [M_x;M_y;M_z] + R_A2*[L(4)/2; 0 ;0] - [U2_x; U2_y; U2_z] - R_U2*[0;0;L(1)/2];
        [U2_x; U2_y; U2_z] + R_U2*[0;0;-L(1)/2] - [L2_x; L2_y; L2_z] - R_L2*[0;0;L(2)/2];
        [L2_x; L2_y; L2_z] + R_L2*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A2*[L(3)/2;0 ;0];  
        ... # Arm 3:  3x3 = 9 constraints
        [M_x;M_y;M_z] + R_A3*[L(4)/2; 0 ;0] - [U3_x; U3_y; U3_z] - R_U3*[0;0;L(1)/2];
        [U3_x; U3_y; U3_z] + R_U3*[0;0;-L(1)/2] - [L3_x; L3_y; L3_z] - R_L3*[0;0;L(2)/2];
        [L3_x; L3_y; L3_z] + R_L3*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A3*[L(3)/2;0 ;0];  
            ... # Fixed variables:  3 constraints
        [M_x;M_y;M_z];
         ... # Fixed lower link spin in own axis
        [0 1 0]*R_L1.'*R_L1_0*[1;0;0];  % L1_y has no component in L1_x_0
        [0 1 0]*R_L2.'*R_L2_0*[1;0;0]; % L2_y has no component in L2_x_0
        [0 1 0]*R_L3.'*R_L3_0*[1;0;0]; % L3_y has no component in L3_x_0
    ];  
    
    % Match R_L symbols to numeric values
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
    
    % Calculate the Jacobian
    J =Jacobian_dNumeric(L(2),...
        rL1_1_1,rL1_1_2,...%R_L1
        rL1_2_1,rL1_2_2,...
        rL1_3_1,rL1_3_2,...
        rL2_1_1,rL2_1_2,...%R_L2
        rL2_2_1,rL2_2_2,...
        rL2_3_1,rL2_3_2,...
        rL3_1_1,rL3_1_2,...%R_L3
        rL3_2_1,rL3_2_2,...
        rL3_3_1,rL3_3_2);        
    % Newton-Raphson update 
        Delta = -J\S;
        qd = qd+Delta; 
        % Update the Rotation matrices by using a constant omega over the
        % time step
                          % Delta theta_L
        deltaR1=expm(tilde(Delta(10:12))); % analytical approximation
        deltaR2=expm(tilde(Delta(19:21)));
        deltaR3=expm(tilde(Delta(28:30)));
        R_L1=deltaR1*R_L1;
        R_L2=deltaR2*R_L2;
        R_L3=deltaR3*R_L3;

        %norm(S) % check the error
    end
    % put results in output variables
    Qd(:,k) = qd;
    R1_out(:,:,k) = R_L1;
    R2_out(:,:,k) = R_L2;
    R3_out(:,:,k) = R_L3;
end
