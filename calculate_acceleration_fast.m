function accelerations = calculate_acceleration_fast(Phi_tt_sym,t_sym,R_L1,R_L2,R_L3,Qvel,thetaA,time_range,L)
    SIZE_Q=42;
    accelerations = zeros(SIZE_Q,length(time_range));
    
    N=length(time_range);
    fprintf('%d time steps. Progess of accelerations: 000.0%%\n',N)
    for timeStep = 1:N
        %fprintf('% d ... %.1f%% complete of accelerations\n',timeStep,100*timeStep/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
        
        Phi_tt=eval(subs(Phi_tt_sym,t_sym,time_range(timeStep))); % Phi_tt should be independent of position
        
        % Match R_L symbols to numeric values
        % R_L1
        R1=R_L1(:,:,timeStep);
        rL1_1_1 = R1(1,1); rL1_1_2 = R1(1,2); rL1_1_3 = R1(1,3);
        rL1_2_1 = R1(2,1); rL1_2_2 = R1(2,2); rL1_2_3 = R1(2,3);
        rL1_3_1 = R1(3,1); rL1_3_2 = R1(3,2); rL1_3_3 = R1(3,3);

        % R_L2
        R2=R_L2(:,:,timeStep);
        rL2_1_1 = R2(1,1); rL2_1_2 = R2(1,2); rL2_1_3 = R2(1,3);
        rL2_2_1 = R2(2,1); rL2_2_2 = R2(2,2); rL2_2_3 = R2(2,3);
        rL2_3_1 = R2(3,1); rL2_3_2 = R2(3,2); rL2_3_3 = R2(3,3);

        % R_L3
        R3=R_L3(:,:,timeStep);
        rL3_1_1 = R3(1,1); rL3_1_2 = R3(1,2); rL3_1_3 = R3(1,3);
        rL3_2_1 = R3(2,1); rL3_2_2 = R3(2,2); rL3_2_3 = R3(2,3);
        rL3_3_1 = R3(3,1); rL3_3_2 = R3(3,2); rL3_3_3 = R3(3,3);
        
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
        
        % Match symbols to numeric velocities
        qvel7=Qvel(7,timeStep); qvel8=Qvel(8,timeStep); qvel9=Qvel(9,timeStep);
        qvel13=Qvel(13,timeStep); qvel14=Qvel(14,timeStep); qvel15=Qvel(15,timeStep);
        qvel19=Qvel(19,timeStep); qvel20=Qvel(20,timeStep); qvel21=Qvel(21,timeStep);
        qvel25=Qvel(25,timeStep); qvel26=Qvel(26,timeStep); qvel27=Qvel(27,timeStep);
        qvel31=Qvel(31,timeStep); qvel32=Qvel(32,timeStep); qvel33=Qvel(33,timeStep);
        qvel37=Qvel(37,timeStep); qvel38=Qvel(38,timeStep); qvel39=Qvel(39,timeStep);
        
        G=GammaNumeric(L(2),L(1),...
            qvel7,qvel8,qvel9,...   %omega_U1
            qvel13,qvel14,qvel15,...%omega_L1
            qvel19,qvel20,qvel21,...%omega_U2
            qvel25,qvel26,qvel27,...%omega_L2
            qvel31,qvel32,qvel33,...%omega_U3
            qvel37,qvel38,qvel39,...%omega_L3
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
        
        accelerations(:,timeStep) = J\(-Phi_tt+G);
    end
end

