function velocities = calculate_velocity_fast( Phi_t_sym,t_sym,R_L1_t,R_L2_t,R_L3_t,thetaA,time_range,L)
    SIZE_Q=42;
    velocities = zeros(SIZE_Q,length(time_range));
    
    % To decrease calculation term further, manually recreate Phi_t
    % But this limits the application, as Phi_t can change easily based on
    % the desired input moments, but the the Jacobian should be the same 
    % for all configurations of the Delta Robot
%     Phi_t=zeros(SIZE_Q,1);
    N=length(time_range);
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
        
        thetaA1=thetaA(1,timeStep);
        thetaA2=thetaA(2,timeStep);
        thetaA3=thetaA(3,timeStep);
        
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