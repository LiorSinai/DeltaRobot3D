function [forces,moments] = calculate_forces(Jacobian,Q_sym,R_sym,thetaA_sym,Q,R_L1_t,R_L2_t,R_L3_t,thetaA,time_range,lambda,indQ,indPhi,r)
    % Use the Lagrange multiplier method to find the forces at a point
    % where the contraint equation in Phi are applied
    forces=zeros(3,length(time_range));
    moments=forces;
    
    % Note unlike in the main code, it is assumed variables are numeric
    % E.g. R_L1 is a numeric matrix
    % unless stated otherwise with _sym
    
    % indQ = indices of the general co-ordinates Q variables. Should be 6
    % consecutives for 3 x,y,z and 3 angles 
    
    % indPhi = indices of the relevant constraints in Phi. Should be two
    % sets of 3 consecutive numbers.
    
    % r = vector from the CoM co-ordinate to the contraint co-ordinate
    
    N=length(time_range);
    fprintf('%d time steps. Progess of forces: 000.0%%\n',N)
    for timeStep = 1:length(time_range)
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
        J=subs(J,thetaA_sym,thetaA(:,timeStep));
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

