function [Q,R1_out,R2_out,R3_out,thetaA_t] = calculate_position_fast(q0,thetaA_0,R_L0,time_range,tolerance,omega,L)
% fast calulation for the position
% does not use the Symbloic toolbox (but that is used to make the Jacobian function)
% This less general than the normal code, as it uses a pre-defined system matrix
% and Jacobian

%% initialise variables
SIZE_Q=42;
Q=zeros(SIZE_Q,length(time_range));
R1_out=zeros(3,3,length(time_range));
R2_out=zeros(3,3,length(time_range));
R3_out=zeros(3,3,length(time_range));
thetaA_t = zeros(3,length(time_range));
singularity=false;
maxIteration=20;

%% extract initial values
R_L1_0=R_L0(:,:,1);
R_L2_0=R_L0(:,:,2);
R_L3_0=R_L0(:,:,3);

thetaU1_y_0=q0(8);  thetaU1_z_0=q0(9);
thetaU2_y_0=q0(20); thetaU2_z_0=q0(21);
thetaU3_y_0=q0(32); thetaU3_z_0=q0(33);

% recalculate R_A matrices
R_A1=rot3D_Rodrigues([0;0;1],thetaU1_z_0); % constant
R_A2=rot3D_Rodrigues([0;0;1],thetaU2_z_0); % constant
R_A3=rot3D_Rodrigues([0;0;1],thetaU3_z_0); % constant

N=length(time_range);
fprintf('%d time steps. Progess of positions: 000.0%%\n',N)
for k = 1:N
    if ~singularity 
        %fprintf('% d ... %.1f%% complete of positions\n',k,100*k/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',k/N*100); %write the new number
    end
    if k == 1
       q=q0;
       R_L1 = R_L1_0; 
       R_L2 = R_L2_0; 
       R_L3 = R_L3_0; 
    end
    % set the error>tolerance on the first loop
    S =ones(SIZE_Q)*100*tolerance;
    newton_iterations = 0;

    while (norm(S) > tolerance)&& (~singularity)
        newton_iterations = newton_iterations +1;
        if newton_iterations>maxIteration
            fprintf('iterions = %d. Probably hit a singularity\n',newton_iterations)
            fprintf('aborting calculations\n')
            % shorten time_range
            time_range=time_range(1:k+1);
            Q=Q(:,1:k+1);
            R1_out=R1_out(:,:,1:k+1);
            R2_out=R2_out(:,:,1:k+1);
            R3_out=R3_out(:,:,1:k+1);
            thetaA_t=thetaA_t(:,1:k+1);
            
            singularity=true;
            break;
        end
        % Match symbols to numeric values
        M_x = q(1);  M_y = q(2);  M_z = q(3);
        U1_x = q(4); U1_y = q(5); U1_z = q(6);
        thetaU1_x=q(7);thetaU1_y=q(8);thetaU1_z=q(9);
        L1_x = q(10); L1_y = q(11); L1_z = q(12);
        thetaL1_x = q(13); thetaL1_y = q(14); thetaL1_z = q(15);
        U2_x = q(16); U2_y = q(17); U2_z = q(18);
        thetaU2_x=q(19);thetaU2_y=q(20);thetaU2_z=q(21);
        L2_x = q(22); L2_y = q(23); L2_z = q(24);
        thetaL2_x = q(25); thetaL2_y = q(26); thetaL2_z = q(27);
        U3_x = q(28); U3_y = q(29); U3_z = q(30);
        thetaU3_x=q(31);thetaU3_y=q(32);thetaU3_z=q(33);
        L3_x = q(34); L3_y = q(35); L3_z = q(36);
        thetaL3_x = q(37); thetaL3_y = q(38); thetaL3_z = q(39);
        P_x = q(40);   P_y = q(41);  P_z = q(42);
        t_sym=time_range(k);

        thetaA1 = thetaU1_x/(-sin(thetaU1_z_0))+ thetaA_0(1);
        thetaA2 = thetaU2_x/(-sin(thetaU2_z_0))+ thetaA_0(2);
        thetaA3 = thetaU3_y;

        thetaA_t(1,k)=thetaA1;
        thetaA_t(2,k)=thetaA2;
        thetaA_t(3,k)=thetaA3;

        R_U1 =R_A1*rot3D_Rodrigues([0;1;0],thetaA1);         
        R_U2 =R_A2*rot3D_Rodrigues([0;1;0],thetaA2); 
        R_U3 =R_A3*rot3D_Rodrigues([0;1;0],thetaA3); 

        S = [...
        ... # Arm 1:  3x3 = 9 constraints
        [M_x;M_y;M_z] + R_A1*[L(4)/2; 0 ;0] - [U1_x; U1_y; U1_z] - R_U1*[0;0;L(1)/2];
        [U1_x; U1_y; U1_z] + R_U1*[0;0;-L(1)/2] - [L1_x; L1_y; L1_z] - R_L1*[0;0;L(2)/2];
        [L1_x; L1_y; L1_z] + R_L1*[0;0;-L(2)/2] - [P_x; P_y; P_z] - R_A1*[L(3)/2;0 ;0];
        ... # upper angles: 3 constraints
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

    J = Jacobian_Numeric(L(2),L(1),...
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

        Delta = -J\S;
        q = q+Delta; % if use q=q-Delta, the code no longer works?
        % Update the Rotation matrices by using a constant omega over the
        % time step
                          % Delta theta_L
%             deltaR1=eye(3)+tilde(Delta(13:15)); % first order approximation
%             deltaR2=eye(3)+tilde(Delta(25:27));
%             deltaR3=eye(3)+tilde(Delta(37:39));
        deltaR1=expm(tilde(Delta(13:15))); % analytical approximation
        deltaR2=expm(tilde(Delta(25:27)));
        deltaR3=expm(tilde(Delta(37:39)));
        R_L1=deltaR1*R_L1;
        R_L2=deltaR2*R_L2;
        R_L3=deltaR3*R_L3;

        %norm(S) % check the error
    end
    % add results to output variables
    Q(:,k) = q;
    R1_out(:,:,k) = R_L1;
    R2_out(:,:,k) = R_L2;
    R3_out(:,:,k) = R_L3;
end

end
