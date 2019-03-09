function [forces,moments] = calculate_forces(Jacobian,Q_sym,R_sym,thetaA_sym,Q,R_L1_t,R_L2_t,R_L3_t,thetaA_t,time_range,lambda,indQ,indPhi,r)
% Calculates forces and moments using the lagrange multiplier method for a
% single point where constraint forces are applied.
% INPUTS
% Jacobian = 42x42 symbolic matrix. Jacobian of the constraint equations
%    Q_sym = 42xN symbolic co-ordinates
%    R_sym = 3x3x6 set of symbolic 3x3 rotation matrices
% thetaA_sym= 3xN symbolic actuator angles
%        Q = 42xN co-ordinate values
% R_L1_t = 3x3xN rotation matrix for lower arm 1
% R_L2_t = 3x3xN rotation matrix for lower arm 2
% R_L3_t = 3x3xN rotation matrix for lower arm 3
% thetaA_t = 3xN actuator angle values
% time_range = 1xN time values
% lamdba = 42xN Lagrange multiplier values
% indQ = indices of the general co-ordinates Q variables. Should be 6
% consecutives for 3 x,y,z and 3 angles
% indPhi = indices of the relevant constraints in Phi. Should be two
% sets of 3 consecutive numbers.
% r = 3xN vector from the CoM co-ordinate to the contraint co-ordinate

% OUTPUTS
% forces =  3xN translational forces
% moments = 3xN moments (rotational forces)

    N=size(thetaA_t,2);
    
    forces=zeros(3,N);
    moments=forces;

    fprintf('%d time steps. Progess of forces: 000.0%%\n',N)
    for timeStep = 1:N
        %fprintf('% d ... %.1f%% complete of forces\n',timeStep,100*timeStep/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete previous number and new line
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
        
        R_L1=R_L1_t(:,:,timeStep);
        R_L2=R_L2_t(:,:,timeStep);
        R_L3=R_L3_t(:,:,timeStep);
        
        
% Much faster to substitute points into the smaller Jacobian
%         J=subs(Jacobian,[R_sym(:,:,4),R_sym(:,:,5),R_sym(:,:,6)],[R_L1,R_L2,R_L3]);
%         J=subs(J,thetaA_sym,thetaA(:,timeStep));
%         J=eval(subs(J,Q_sym,Q(:,timeStep)));
        J=Jacobian(indPhi,indQ);
        J=subs(J,[R_sym(:,:,4),R_sym(:,:,5),R_sym(:,:,6)],[R_L1,R_L2,R_L3]);
        J=subs(J,thetaA_sym,thetaA_t(:,timeStep));
        J=eval(subs(J,Q_sym(indQ),Q(indQ,timeStep)));
        
        % equation 6.43 on pg 71 of the notes
        lambda_i=lambda(indPhi,timeStep);
%         Phi_ri    =J(indPhi,indQ(1:3));
%         Phi_thetai=J(indPhi,indQ(4:6));
        Phi_ri    =J(:,1:3);
        Phi_thetai=J(:,4:6);   
        forces(:,timeStep)   =-Phi_ri.'*lambda_i;
        moments(:,timeStep)=-Phi_thetai.'*lambda_i-tilde(r(:,timeStep))*forces(:,timeStep);
    end
end

