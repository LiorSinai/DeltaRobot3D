function velocities = calculate_velocity( Velocity,Q_sym,R_sym,thetaA_sym,t_sym,Q,R_L1,R_L2,R_L3,thetaA,time_range)
% Calculate velocities for the co-ordinates given known positions.   

% INPUTS
% Velocity = 42x1 symbolic vector for velocities based on Phi.
%  Q_sym = 42xN symbolic co-ordinates
%  R_sym = (3x3)x6 set of symbolic 3x3 rotation matrices
% thetaA_sym= 3xN symbolic actuator angles
%  t_sym = 1x1 time symbol
%   R_L1 = 3x3xN rotation matrix for lower arm 1
%   R_L2 = 3x3xN rotation matrix for lower arm 2
%   R_L3 = 3x3xN rotation matrix for lower arm 3
%   R_L0 = (3x3)x3 set of 3x3 rotation matrices. Initial values (guess)
% thetaA = 3xN actuator angle values
% time_range = 1xN time values

% OUTPUTS
% velocities =42xN velocities for each co-ordinate at each time step

    N=size(thetaA,2);
    velocities = zeros(length(Q_sym),N);
    
    fprintf('%d time steps. Progess of velocities: 000.0%%\n',N)
    for timeStep = 1:N
        %fprintf('% d ... %.1f%% complete of velocities\n',step,100*step/length(time_range))
        fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
        fprintf('%05.1f%%\n',timeStep/N*100); %write the new number
        
        temp=subs(Velocity,t_sym,time_range(timeStep));
        temp=subs(temp,[R_sym(:,:,4),R_sym(:,:,5),R_sym(:,:,6)],[R_L1(:,:,timeStep),R_L2(:,:,timeStep),R_L3(:,:,timeStep)]);
        temp=subs(temp,thetaA_sym,thetaA(:,timeStep));
        velocities(:,timeStep) = eval(subs(temp,Q_sym,Q(:,timeStep)));
        
    end
end