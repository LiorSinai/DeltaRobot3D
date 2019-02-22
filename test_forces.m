% 20 February 2019
% Test forces
% Should run prokect_week_5.m and project_week_5.force first to get:
%   all symbols
%   all rotation matrices
%   Q_AC
%   the force structure

%% check errors
errorForcesResidue=zeros(9,length(time_range));
errorMomentsResidue=zeros(9,length(time_range));

for timeStep = 1:length(time_range)
    % map values
    forcesA1=wrenches(1).forces(:,timeStep);
    forcesB1=wrenches(2).forces(:,timeStep);
    forcesA2=wrenches(3).forces(:,timeStep);
    forcesB2=wrenches(4).forces(:,timeStep);
    forcesA3=wrenches(5).forces(:,timeStep);
    forcesB3=wrenches(6).forces(:,timeStep);
    
    % force balance in the inertial frame 0
    % m_u*rddot =FA+FB+[0;0;-m_u*g]
    errorForcesResidue(1:3,timeStep)=m_u*Qacc(4:6,timeStep)+...
                               -forcesA1-forcesB1-[0;0;-m_u*g];
   errorForcesResidue(4:6,timeStep)=m_u*Qacc(16:18,timeStep)+...
                               -forcesA2-forcesB2-[0;0;-m_u*g];
   errorForcesResidue(7:9,timeStep)=m_u*Qacc(28:30,timeStep)+...
                               -forcesA3-forcesB3-[0;0;-m_u*g];
   
   % moment balance in the body fixed frame C
   % I_c*alphaC+tilde(omegaC)*I_c*omegaC = tilde(r^U_A)*F^A_U+tilde(r^U_B)*F^A_B
   %                                       +[0;M_A;0]
   R=eval(subs(R_U1,thetaA1,thetaA_t(1,timeStep)));
   alphaC=R'*Qacc(7:9,timeStep);
   omegaC=R'*Qvel(7:9,timeStep); 
   momentsC=wrenches(1).momentsC(:,timeStep);
   errorMomentsResidue(1:3,timeStep)=I_C(1:3).*R'*Qacc(7:9,timeStep)+...
                                +tilde(omegaC)*diag(I_C(1:3))*omegaC+...
                                -tilde([0; 0 ;+L(1)/2])*R'*forcesA1+...
                                -tilde([0; 0 ;-L(1)/2])*R'*forcesB1+...
                                -momentsC;
   % components of the error
   M1(:,timeStep)=I_C(1:3).*alphaC;
   M2(:,timeStep)=tilde(omegaC)*diag(I_C(1:3))*omegaC;
   M3(:,timeStep)=tilde([0; 0 ;+L(1)/2])*R'*forcesA1;
   M4(:,timeStep)=tilde([0; 0 ;-L(1)/2])*R'*forcesB1;
   M5(:,timeStep)=momentsC;
   
   R=eval(subs(R_U2,thetaA2,thetaA_t(2,timeStep)));
   omegaC=R'*Qvel(19:21,timeStep); 
   alphaC=R'*Qacc(19:21,timeStep);
   momentsC=wrenches(3).momentsC(:,timeStep);
   errorMomentsResidue(4:6,timeStep)=I_C(1:3).*R'*alphaC+...
                                +tilde(omegaC)*diag(I_C(1:3))*omegaC+...
                                -tilde([0; 0 ;+L(1)/2])*R'*forcesA2+...
                                -tilde([0; 0 ;-L(1)/2])*R'*forcesB2+...
                                -momentsC;

   R=eval(subs(R_U3,thetaA3,thetaA_t(3,timeStep)));
   omegaC=R'*Qvel(31:33,timeStep); 
   alphaC=R'*Qacc(31:33,timeStep);
   momentsC=wrenches(5).momentsC(:,timeStep);
   errorMomentsResidue(7:9,timeStep)=I_C(1:3).*alphaC+...
                                +tilde(omegaC)*diag(I_C(1:3))*omegaC+...
                                -tilde([0; 0 ;+L(1)/2])*R'*forcesA3+...
                                -tilde([0; 0 ;-L(1)/2])*R'*forcesB3+...
                                -momentsC;
                            
end

%% plot components
% Note: for twist, there will be an applied moment about both the body frame
% x-axis and y-axis for arms 2 and 3.
% This does not seem to be an error, as the asymmetry of the sitatuion
% tends to rotate the delta robot about the x axis, and therefore
% an applied moment has to counter-act this
N=length(time_range);
figure;plot(1:N,M1,'-xb',1:N,M2,'-xr',1:N,M3,'-xg',1:N,M4,'-xm',1:N,M5,'-xk');
grid;legend;

%% plot error residues
figure; 
plot(1:length(time_range),errorForcesResidue,'-x'); legend; grid
title('Error residues for forces')
figure; 
plot(1:length(time_range),errorMomentsResidue,'-x'); legend; grid
title('Error residues for moments')

%save(sprintf('TestRun_%s',timestamp)); % save all variables

