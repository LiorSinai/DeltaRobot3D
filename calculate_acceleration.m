function accelerations = calculate_acceleration(Jacobian,Phi_tt_sym,Q_sym,R_sym,thetaA_sym,t_sym,Q,R_L1_t,R_L2_t,R_L3_t,Qvel,thetaA_t,time_range,L)
% Caluculates accelerations given all velocities and positions. Mostly
% numeric functions. Only Phi_tt uses symbolic toolbox functions.
% INPUTS
% Jacobian = 42x42 symbolic Jacobian of the system wrt 42 co-ordinates
% Phi_tt_sym = 42x1 symbolic function, second time derivative of the constraint equations
%  Q_sym = 42xN symbolic co-ordinates
%  R_sym = (3x3)x6 set of symbolic 3x3 rotation matrices
% thetaA_sym= 3xN symbolic actuator angles
%  t_sym = 1x1 time symbol
%      Q = 42xN co-ordinate positions
% R_L1_t = (3x3)xN rotation matrix for lower arm 1
% R_L2_t = (3x3)xN rotation matrix for lower arm 2
% R_L3_t = (3x3)xN rotation matrix for lower arm 3
% Qvel   = 42xN velocites
% thetaA_t = 3xN actuator angles
% time_range = 1xN time values
% L=[L_upper L_lower L_endEffector L_base] ... lengths [m]

% OUTPUTS
% accelerations = 42xN accelerations 

    N=size(thetaA_t,2);
    accelerations = zeros(length(Q_sym),N);
    
   R_L1_0=R_L1_t(:,:,1);
   R_L2_0=R_L2_t(:,:,1);
   R_L3_0=R_L3_t(:,:,1);
    
    fprintf('%d time steps. Progess of accelerations: 000.0%%\n',N)
    for timeStep = 1:N
        %fprintf('% d ... %.1f%% complete of accelerations\n',timeStep,100*timeStep/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
        
        Vt=eval(subs(Phi_tt_sym,t_sym,time_range(timeStep)));
        
        R_U1=eval(subs(R_sym(:,:,1),thetaA_sym(1),thetaA_t(1,timeStep)));
        R_U2=eval(subs(R_sym(:,:,2),thetaA_sym(2),thetaA_t(2,timeStep)));
        R_U3=eval(subs(R_sym(:,:,3),thetaA_sym(3),thetaA_t(3,timeStep)));
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
        J=subs(J,thetaA_sym,thetaA_t(:,timeStep));
        J=eval(subs(J,Q_sym,Q(:,timeStep)));
        
        accelerations(:,timeStep) = J\(Vt-Gamma);
        %accelerations(:,time+1) = eval(subs(subs(Acceleration(:,time+1),Q,Phi(:,time+1)),[R_L1,R_L2,R_L3],[R1(:,:,time+1),R2(:,:,time+1),R3(:,:,time+1)]));
        %velocities(:,time+1) = eval(subs(subs(Velocity,[R_L1,R_L2,R_L3],[R1(:,:,time+1),R2(:,:,time+1),R3(:,:,time+1)]),Q,Phi(:,time+1)));
    end
end

