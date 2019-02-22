function accelerations = calculate_acceleration(Jacobian,Vt_sym,Q_sym,R_sym,thetaA_sym,t_sym,Q,R_L1_t,R_L2_t,R_L3_t,Qvel,thetaA,time_range,L)
    accelerations = zeros(length(Q_sym),length(time_range));
    
    % Note unlike in the main code, it is assumed variables are numeric
    % E.g. R_L1 is a numeric matrix
    % unless stated otherwise with _sym
   R_L1_0=R_L1_t(:,:,1);
   R_L2_0=R_L2_t(:,:,1);
   R_L3_0=R_L3_t(:,:,1);
    
    N=length(time_range);
    fprintf('%d time steps. Progess of accelerations: 000.0%%\n',N)
    for timeStep = 1:length(time_range)
        %fprintf('% d ... %.1f%% complete of accelerations\n',timeStep,100*timeStep/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
        
        Vt=eval(subs(Vt_sym,t_sym,time_range(timeStep)));
        
        R_U1=eval(subs(R_sym(:,:,1),thetaA_sym(1),thetaA(1,timeStep)));
        R_U2=eval(subs(R_sym(:,:,2),thetaA_sym(2),thetaA(2,timeStep)));
        R_U3=eval(subs(R_sym(:,:,3),thetaA_sym(3),thetaA(3,timeStep)));
        R_L1=R_L1_t(:,:,timeStep);
        R_L2=R_L2_t(:,:,timeStep);
        R_L3=R_L3_t(:,:,timeStep);
        
        Gamma_1 = zeros(size(Jacobian)); % reset Gamma_1 to zero on every loop
        Gamma_1(1:3,7:9)    =-tilde(Qvel(7:9,timeStep))*R_U1*tilde([0;0;+L(1)/2]).'*R_U1.';
        Gamma_1(4:6,7:9)    =+tilde(Qvel(7:9,timeStep))*R_U1*tilde([0;0;-L(1)/2]).'*R_U1.';
        Gamma_1(4:6,13:15)  =-tilde(Qvel(13:15,timeStep))*R_L1*tilde([0;0;+L(2)/2]).'*R_L1.';
        Gamma_1(7:9,13:15)  =+tilde(Qvel(13:15,timeStep))*R_L1*tilde([0;0;-L(2)/2]).'*R_L1.';

        Gamma_1(13:15,19:21) =-tilde(Qvel(19:21,timeStep))*R_U2*tilde([0;0;+L(1)/2]).'*R_U2.';
        Gamma_1(16:18,19:21) =+tilde(Qvel(19:21,timeStep))*R_U2*tilde([0;0;-L(1)/2]).'*R_U2.';
        Gamma_1(16:18,25:27) =-tilde(Qvel(25:27,timeStep))*R_L2*tilde([0;0;+L(2)/2]).'*R_L2.';
        Gamma_1(19:21,25:27) =+tilde(Qvel(25:27,timeStep))*R_L2*tilde([0;0;-L(2)/2]).'*R_L2.';

        Gamma_1(25:27,31:33) =-tilde(Qvel(31:33,timeStep))*R_U3*tilde([0;0;+L(1)/2]).'*R_U3.';
        Gamma_1(28:30,31:33) =+tilde(Qvel(31:33,timeStep))*R_U3*tilde([0;0;-L(1)/2]).'*R_U3.';
        Gamma_1(28:30,37:39) =-tilde(Qvel(37:39,timeStep))*R_L3*tilde([0;0;+L(2)/2]).'*R_L3.';
        Gamma_1(31:33,37:39) =+tilde(Qvel(37:39,timeStep))*R_L3*tilde([0;0;-L(2)/2]).'*R_L3.';

        %lower link Jacobian
        Gamma_1(40,13:15) = [0 1 0]*(tilde(Qvel(13:15,timeStep))*R_L1).'*tilde(R_L1_0*[1;0;0]);
        Gamma_1(41,25:27) = [0 1 0]*(tilde(Qvel(25:27,timeStep))*R_L2).'*tilde(R_L2_0*[1;0;0]);
        Gamma_1(42,37:39) = [0 1 0]*(tilde(Qvel(37:39,timeStep))*R_L3).'*tilde(R_L3_0*[1;0;0]);

        Gamma = Gamma_1*Qvel(:,timeStep);
        
        J=subs(Jacobian,[R_sym(:,:,4),R_sym(:,:,5),R_sym(:,:,6)],[R_L1,R_L2,R_L3]);
        J=subs(J,thetaA_sym,thetaA(:,timeStep));
        J=eval(subs(J,Q_sym,Q(:,timeStep)));
        
        accelerations(:,timeStep) = J\(Vt-Gamma);
        %accelerations(:,time+1) = eval(subs(subs(Acceleration(:,time+1),Q,Phi(:,time+1)),[R_L1,R_L2,R_L3],[R1(:,:,time+1),R2(:,:,time+1),R3(:,:,time+1)]));
        %velocities(:,time+1) = eval(subs(subs(Velocity,[R_L1,R_L2,R_L3],[R1(:,:,time+1),R2(:,:,time+1),R3(:,:,time+1)]),Q,Phi(:,time+1)));
    end
end

