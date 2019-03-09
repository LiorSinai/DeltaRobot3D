function velocities = calculate_velocity_fast( Phi_t_sym,t_sym,R_L1_t,R_L2_t,R_L3_t,thetaA_t,time_range,L)
% Calculate velocities for the co-ordinates given known positions.  This
% version reduces reliance on the symbolic toolbox to a minimum, which 
%improves the performance.

% INPUTS
%  Phi_t_sym = 42x1 symbolic derivative of the constrain equation vector Phi
%  t_sym = 1x1 time symbol
% thetaA_sym= 3xN symbolic actuator angles
%   R_L1_t = 3x3xN rotation matrix for lower arm 1
%   R_L2_t = 3x3xN rotation matrix for lower arm 2
%   R_L3_t = 3x3xN rotation matrix for lower arm 3
% thetaA_t = 3xN actuator angle values
% time_range = 1xN time values
% L=[L_upper L_lower L_endEffector L_base] ... lengths [m]

% OUTPUTS
% velocities =42xN velocities for each co-ordinate at each time step    

    SIZE_Q=42;
    N=size(thetaA_t,2);
    velocities = zeros(SIZE_Q,N);
    
    % Manually recreate Phi_t if not using the symbolic toolbox:
%     Phi_t=zeros(SIZE_Q,1);

    fprintf('%d time steps. Progess of velocities: 000.0%%\n',N)
    for timeStep = 1:N
        %fprintf('% d ... %.1f%% complete of velocities\n',step,100*step/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
        
        Phi_t=eval(subs(Phi_t_sym,t_sym,time_range(timeStep)));
        % Phi_t should not be dependent on the co-ordinates
        
        % Match R_L symbols to numeric values
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
        
        thetaA1=thetaA_t(1,timeStep);
        thetaA2=thetaA_t(2,timeStep);
        thetaA3=thetaA_t(3,timeStep);
        
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

        velocities(:,timeStep) =-J\Phi_t;
        
    end
end