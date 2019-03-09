function [Q,R1_out,R2_out,R3_out,thetaA_t] = calculate_position( System,Jacobian,Q_sym,R_sym,thetaA_sym,t_sym,q0,thetaA_0,R_L0,time_range,tolerance )
%Calcualtes the position co-ordiantes of a 3D Delta robot. Relies heavily
%on the symbolic toolbox and therefore is very slow.

% INPUTS
% System = 42x1 symbolic system of constraint equations
% Jacobian = 42x42 symbolic Jacobian of the system wrt 42 co-ordinates
%  Q_sym = 42xN symbolic co-ordinates
%  R_sym = 3x3x6 set of symbolic 3x3 rotation matrices
% thetaA_sym= 3xN symbolic actuator angles
%  t_sym = 1x1 symbol for time
%     q0 = 42x1 initial position (guess)
%thetaA_0 = 3x1 initial actuator angles (guess) 
%   R_L0 = (3x3)x3 set of 3x3 rotation matrices. Initial values (guess)
% time_range = 1xN time values. The first value should be 0 to confirm the
% intial values/guesses
% tolerance = desired accuracy for the Newton-Raphson iterations

% OUTPUTS
% Note: if a singularity is hit, N is changed to N=singularity time step
% Q = 42xN co-ordinate values
% R1_out= 3x3xN rotation matrices set for lower arm 1
% R2_out= 3x3xN rotation matrices set for lower arm 2
% R3_out= 3x3xN rotation matrices set for lower arm 3
% thetaA_t = 3xN actuator angles

Q=zeros(length(Q_sym),length(time_range));
R1_out=zeros(3,3,length(time_range));
R2_out=zeros(3,3,length(time_range));
R3_out=zeros(3,3,length(time_range));
thetaA_t = zeros(3,length(time_range));
singularity=false;
maxIteration=20;

R_U1=R_sym(:,:,1);
R_U2=R_sym(:,:,2);
R_U3=R_sym(:,:,3);
R_L1=R_sym(:,:,4);
R_L2=R_sym(:,:,5);
R_L3=R_sym(:,:,6);

R_L1_0=R_L0(:,:,1);
R_L2_0=R_L0(:,:,2);
R_L3_0=R_L0(:,:,3);

N=length(time_range);
fprintf('%d time steps. Progess of positions: 000.0%%\n',N)
    for k = 1:N
        if ~singularity 
            %fprintf('% d ... %.1f%% complete of positions\n',k,100*k/length(time_range))
            fprintf('\b\b\b\b\b\b\b')  %delete new line, % sign, previous number
            fprintf('%05.1f%%\n',k/N*100); %write the new number
        end
        if k == 1
           q=q0;
           R1 = R_L1_0; 
           R2 = R_L2_0; 
           R3 = R_L3_0; 
        end
        %Phi(:,time) = q;
        % set the error>tolerance on the first loop
        S =ones(size(Q_sym))*100*tolerance;
        newton_iterations = 0;
        
        while (norm(S) > tolerance)&& (~singularity)
            newton_iterations = newton_iterations +1;
            if newton_iterations>maxIteration
                warning('iterions = %d. Probably hit a singularity\n',newton_iterations)
                warning('aborting calculations\n')
                % shorten time_range
                time_range=time_range(1:k+1);
                singularity=true;
                break;
            end
            thetaA_t(1,k) = q(7)/R_U1(1,2) + thetaA_0(1);
            thetaA_t(2,k) = q(19)/R_U2(1,2)+ thetaA_0(2);
            thetaA_t(3,k) = q(32);
            
            S =  subs(System,[Q_sym;t_sym],[q;time_range(k)]);
            J =  subs(Jacobian,[Q_sym;t_sym],[q;time_range(k)]); 
            
            S = subs(S,thetaA_sym,thetaA_t(:,k));
            J = subs(J,thetaA_sym,thetaA_t(:,k));
                                  
            S=eval(subs(S,[R_L1,R_L2,R_L3],[R1,R2,R3]));
            J=eval(subs(J,[R_L1,R_L2,R_L3],[R1,R2,R3]));
          
            Delta = -J\S;
            q = q+Delta; % if use q=q-Delta, the code no longer works?
            % Update the Rotation matrices by using a constant omega over the
            % time step
                             % Delta theta_L
%             deltaR1=eye(3)+tilde(Delta(13:15)); % first order approximation
%             deltaR2=eye(3)+tilde(Delta(25:27));
%             deltaR3=eye(3)+tilde(Delta(37:39));
            deltaR1=expm(tilde(Delta(13:15))); % analytical approximation
            deltaR2=expm(tilde(Delta(25:27)));
            deltaR3=expm(tilde(Delta(37:39)));
            R1=deltaR1*R1;
            R2=deltaR2*R2;
            R3=deltaR3*R3;
            
            %norm(S)
        end
        Q(:,k) = q;
        R1_out(:,:,k) = R1;
        R2_out(:,:,k) = R2;
        R3_out(:,:,k) = R3;
    end
end
