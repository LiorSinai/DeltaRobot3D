function [forces,moments] = calculate_forces_dependent(R_L1_t,R_L2_t,R_L3_t,thetaA_t,L,lambda,indQd,indPhi,r)
% Use the Lagrange multiplier method to find the forces at a point
% where the contraint equation in Phi are applied
    
% INPUTS
% R_L1_t = 3x3xN rotation matrix for lower arm 1
% R_L2_t = 3x3xN rotation matrix for lower arm 2
% R_L3_t = 3x3xN rotation matrix for lower arm 3
% thetaA_t = 3xN actuator angle values
% L=[L_upper L_lower L_endEffector L_base] ... lengths
% lamdba = 42xN Lagrange multiplier values   
% indQd = 6x1 indices of the dependent co-ordinates. Should be 6
% consecutives for 3 x,y,z and 3 angles, or at least 1 NaN for the
% independent co-ordinates.
% indPhi = 6x1 indices of the relevant constraints in Phi_Without_driving. 
%Should be two sets of 3 consecutive numbers, or at least 1 NaN for the 
%indepedent co-ordinates, because there are no constraint equations in Phi 
%for the indepedent co-ordinates.
% r = 3xN vector from the CoM co-ordinate to the contraint co-ordinate

% OUTPUTS
% forces =  3xN translational forces in the inertial frame
% moments = 3xN moments (rotational forces) in the inertial frame
    
    N=size(thetaA_t,2);
    
    forces=zeros(3,N);
    moments=forces;    
    
    fprintf('%d time steps. Progess of forces: 000.0%%\n',N)
    for timeStep = 1:N
        %fprintf('% d ... %.1f%% complete of forces\n',timeStep,100*timeStep/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete previous number and new line
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
              
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

        J = Jacobian_dNumeric(L(2),...
        rL1_1_1,rL1_1_2,...%R_L1
        rL1_2_1,rL1_2_2,...
        rL1_3_1,rL1_3_2,...
        rL2_1_1,rL2_1_2,...%R_L2
        rL2_2_1,rL2_2_2,...
        rL2_3_1,rL2_3_2,...
        rL3_1_1,rL3_1_2,...%R_L3
        rL3_2_1,rL3_2_2,...
        rL3_3_1,rL3_3_2);
    
        if sum(isnan(indPhi)>0)
            % cannot find forces for the independent co-ordinates
            % therefore do not do that calculation
            indValid=~isnan(indPhi);
            lambda_i=lambda(indPhi(indValid),timeStep);
            Phi_ri    =J(indPhi(indValid),indQd(indValid));
            forces(:,timeStep)   =-Phi_ri.'*lambda_i;
        else      
            % equation 6.43 on pg 71 of the notes
            lambda_i=lambda(indPhi,timeStep);
            Phi_ri    =J(indPhi,indQd(1:3));
            Phi_thetai=J(indPhi,indQd(4:6));
            forces(:,timeStep)   =-Phi_ri.'*lambda_i;
            moments(:,timeStep)=-Phi_thetai.'*lambda_i-tilde(r(:,timeStep))*forces(:,timeStep);
        end
    end
end

