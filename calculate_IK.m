%% Inverse Kinematics
%Lior Sinai, 2019-02-24

function [qa, qaVel, qaAcc]=calculate_IK(L,qa0,P,Pvel,Pacc,tolerance)
% settings
singularity=false;
maxIteration=20;

% extract variables
L_u=L(1);
L_l=L(2);
L_e=L(3);
L_b=L(4);
P_x=P(1);P_y=P(2);P_z=P(3);
Pvel_x=Pvel(1);Pvel_y=Pvel(2); Pvel_z=Pvel(3);
Pacc_x=Pacc(1);Pacc_y=Pacc(2); Pacc_z=Pacc(3);

qa = zeros(3,1);
qaAcc = zeros(3,1);
thetaU1_z_0=-2*pi/3;
thetaU2_z_0=+2*pi/3;

Phi_IK =ones(3)*100*tolerance;
iter = 0;
%% Position
while (norm(Phi_IK) > tolerance)&& (~singularity)
    iter = iter +1;
    if iter>maxIteration
        warning('iterions = %d. Probably hit a singularity\n',iter)
        warning('aborting calculations\n')
        % shorten time_range
        singularity=true;
        break;
    end
    if iter==1
        qa=qa0; % set initial value
    end
    % Map variables
    thetaA1=qa(1);
    thetaA2=qa(2);
    thetaA3=qa(3);
    
    Phi_IK=[
        ((0.5*(L_b-L_e)-L_u*sin(thetaA1))*cos(thetaU1_z_0)-P_x)^2+...
        ((0.5*(L_b-L_e)-L_u*sin(thetaA1))*sin(thetaU1_z_0)-P_y)^2+...
        (-L_u*cos(thetaA1)-P_z)^2-L_l^2;
        ((0.5*(L_b-L_e)-L_u*sin(thetaA2))*cos(thetaU2_z_0)-P_x)^2+...
        ((0.5*(L_b-L_e)-L_u*sin(thetaA2))*sin(thetaU2_z_0)-P_y)^2+...
        (-L_u*cos(thetaA2)-P_z)^2-L_l^2;
        (0.5*(L_b-L_e)-L_u*sin(thetaA3)-P_x)^2+(-P_y)^2+...
        (-L_u*cos(thetaA3)-P_z)^2-L_l^2
        ];
    Jacobian_qa = Jacobian_qaNumeric(L_b,L_e,L_u,P_x,P_y,P_z,thetaA1,thetaA2,thetaA3);
    Delta = -Jacobian_qa\Phi_IK;
    qa = qa+Delta; 
end
% extract final variables
thetaA1=qa(1);
thetaA2=qa(2);
thetaA3=qa(3);

%% Velocity
Jacobian_P=Jacobian_PNumeric(L_b,L_e,L_u,...
    P_x,P_y,P_z,...
    thetaA1,thetaA2,thetaA3);
qaVel=-Jacobian_qa\Jacobian_P*Pvel;

%% Acceleration
omegaA1=qaVel(1);
omegaA2=qaVel(2);
omegaA3=qaVel(3);
gamma_IK = Gamma_IKNumeric(L_b,L_e,L_u,...
    P_x,P_y,P_z,...
    Pacc_x,Pacc_y,Pacc_z,...
    Pvel_x,Pvel_y,Pvel_z,...
    omegaA1,omegaA2,omegaA3,...
    thetaA1,thetaA2,thetaA3);
qaAcc=Jacobian_qa\gamma_IK;    
  