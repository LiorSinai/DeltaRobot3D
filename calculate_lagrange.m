function [lambda, Q_AC] = calculate_lagrange(Jacobian,Q_sym,R_sym,thetaA_sym,Q,Qacc,R_L1_t,R_L2_t,R_L3_t,Qvel,thetaA,time_range,M_C,I_C,Q_A)
    lambda=zeros(size(Q));
    Q_AC=zeros(size(Q));
    
    % Note unlike in the main code, it is assumed variables are numeric
    % E.g. R_L1 is a numeric matrix
    % unless stated otherwise with _sym
   
   I_upper_xx=I_C(1);
   I_upper_yy=I_C(2);
   I_upper_zz=I_C(3);
   I_lower_xx=I_C(4);
   I_lower_yy=I_C(5);
   I_lower_zz=I_C(6);
      
    N=length(time_range);
    fprintf('%d time steps. Progess of lambdas: 000.0%%\n',N)
    for timeStep = 1:length(time_range)
        %fprintf('% d ... %.1f%% complete of lambdas\n',timeStep,100*timeStep/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
        
        R_U1=eval(subs(R_sym(:,:,1),thetaA_sym(1),thetaA(1,timeStep)));
        R_U2=eval(subs(R_sym(:,:,2),thetaA_sym(2),thetaA(2,timeStep)));
        R_U3=eval(subs(R_sym(:,:,3),thetaA_sym(3),thetaA(3,timeStep)));
        R_L1=R_L1_t(:,:,timeStep);
        R_L2=R_L2_t(:,:,timeStep);
        R_L3=R_L3_t(:,:,timeStep);

        M_0=M_C;  % reset mass (moment of inertia) matrix
        Q_AC(:,timeStep)=Q_A; % reset applied forces part

        % Calculate moment inertias in the inertia frame
        I_U1=R_U1*diag([I_upper_xx,I_upper_yy, I_upper_zz])*R_U1.';
        I_U2=R_U2*diag([I_upper_xx,I_upper_yy, I_upper_zz])*R_U2.';
        I_U3=R_U3*diag([I_upper_xx,I_upper_yy, I_upper_zz])*R_U3.';

        I_L1=R_L1*diag([I_lower_xx,I_lower_yy, I_lower_zz])*R_L1.';
        I_L2=R_L2*diag([I_lower_xx,I_lower_yy, I_lower_zz])*R_L2.';
        I_L3=R_L3*diag([I_lower_xx,I_lower_yy, I_lower_zz])*R_L3.';

        % correct moment of inertias
        M_0(7:9,7:9)    =I_U1;
        M_0(13:15,13:15)=I_L1;    
        M_0(19:21,19:21)=I_U2;  
        M_0(25:27,25:27)=I_L2;  
        M_0(31:33,31:33)=I_U3;
        M_0(37:39,37:39)=I_L3;

        % Calculate gyroscopic forces in the inertia frame
        C_U1=tilde(Qvel(7:9,timeStep))*  I_U1*Qvel(7:9,timeStep);
        C_U2=tilde(Qvel(19:21,timeStep))*I_U2*Qvel(19:21,timeStep);
        C_U3=tilde(Qvel(31:33,timeStep))*I_U3*Qvel(31:33,timeStep);

        C_L1=tilde(Qvel(13:15,timeStep))*I_L1*Qvel(13:15,timeStep);
        C_L2=tilde(Qvel(25:27,timeStep))*I_L2*Qvel(25:27,timeStep);
        C_L3=tilde(Qvel(37:39,timeStep))*I_L3*Qvel(37:39,timeStep);

        % Q_A_C = Q_A-C*qdot
        Q_AC(7:9,timeStep)  =Q_A(7:9)  -C_U1;
        Q_AC(13:15,timeStep)=Q_A(13:15)-C_L1;
        Q_AC(19:21,timeStep)=Q_A(19:21)-C_U2;
        Q_AC(25:27,timeStep)=Q_A(25:27)-C_L2;
        Q_AC(31:33,timeStep)=Q_A(31:33)-C_U3;
        Q_AC(37:39,timeStep)=Q_A(37:39)-C_L3;

        J=subs(Jacobian,[R_sym(:,:,4),R_sym(:,:,5),R_sym(:,:,6)],[R_L1,R_L2,R_L3]);
        J=subs(J,thetaA_sym,thetaA(:,timeStep));
        J=eval(subs(J,Q_sym,Q(:,timeStep)));

        lambda(:,timeStep)=J.'\(Q_AC(:,timeStep)-M_0*Qacc(:,timeStep));
    end
end

