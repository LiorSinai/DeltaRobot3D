function plot_velocity( velocities,time_range )
% plots the components of each centre of mass' velocity seperately wrt
% time. Plots are in a 3x2 subfigure and a seperate figure for the end 
% effector

% INPUTS 
% velocities = 42xN velocity values
% time_range = 1xN time values

% Q=[M_x;  M_y;  M_z; 
%     % arm 1
%    U1_x; U1_y; U1_z; thetaU1_x; thetaU1_y; thetaU1_z;
%    L1_x; L1_y; L1_z; thetaL1_x; thetaL1_y; thetaL1_z;
%     % arm 2
%    U2_x; U2_y; U2_z; thetaU2_x; thetaU2_y; thetaU2_z;
%    L2_x; L2_y; L2_z; thetaL2_x; thetaL2_y; thetaL2_z;
%     % arm 3
%    U3_x; U3_y; U3_z; thetaU3_x; thetaU3_y; thetaU3_z;
%    L3_x; L3_y; L3_z; thetaL3_x; thetaL3_y; thetaL3_z;
%    P_x;  P_y; P_z
%   ];

    % centers of mass arm1
    U1vx = velocities(4,:);
    U1vy = velocities(5,:);
    U1vz = velocities(6,:);
    omegaU1=velocities(7:9,:);
    
    L1vx = velocities(10,:);
    L1vy = velocities(11,:);
    L1vz = velocities(12,:);
    omegaL1=velocities(13:15,:);
    
    % second center of mass
    U2vx = velocities(16,:);
    U2vy = velocities(17,:);
    U2vz = velocities(18,:);
    omegaU2=velocities(19:21,:);
    
    L2vx = velocities(22,:);
    L2vy = velocities(23,:);
    L2vz = velocities(24,:);
    omegaL2=velocities(25:27,:);
    
    % third center of mass
    U3vx = velocities(28,:);
    U3vy = velocities(29,:);
    U3vz = velocities(30,:);
    omegaU3=velocities(31:33,:);

    L3vx = velocities(34,:);
    L3vy = velocities(35,:);
    L3vz = velocities(36,:);
    omegaL3=velocities(37:39,:);
       
    %end-effector
    
    Pxv = velocities(40,:);
    Pyv = velocities(41,:);
    Pzv = velocities(42,:);
    
   
    figure('Name','Velocity');

    %plot velocities of 2 centers of mass arm1 
    subplot(3,2,1)
    hold on
    plot(time_range,U1vx,'-')
    plot(time_range,U1vy,'-')
    plot(time_range,U1vz,'--')
    plot(time_range,L1vx,'--')
    plot(time_range,L1vy,'--')
    plot(time_range,L1vz,'--')
    legend('U1v_x','U1v_y','U1v_z','L1v_x','L1v_y','L1v_z')
    title('Centre of masses of arm 1: linear velocities');
    grid;
    subplot(3,2,2)
    hold on
    plot(time_range,omegaU1,'-')
    plot(time_range,omegaL1,'--')
    title('Centre of masses of arm 1: angular velocites');
    legend('\omega_{U1,x}','\omega_{U1,y}','\omega_{U1,z}',...
           '\omega_{L1,x}','\omega_{L1,y}','\omega_{L1,z}')
    grid;
    
    %plot velocities of 2 centers of mass arm2 
    subplot(3,2,3)
    hold on
    plot(time_range,U2vx,'-')
    plot(time_range,U2vy,'-')
    plot(time_range,U2vz,'--')
    plot(time_range,L2vx,'--')
    plot(time_range,L2vy,'--')
    plot(time_range,L2vz,'--')
    legend('U2v_x','U2v_y','U2v_z','L2v_x','L2v_y','L2v_z')
    title('Centre of masses of arm 2: linear velocities');
    grid;
    subplot(3,2,4)
    hold on
    plot(time_range,omegaU2,'-')
    plot(time_range,omegaL2,'--')
    title('Centre of masses of arm 2: angular velocites');
    legend('\omega_{U2,x}','\omega_{U2,y}','\omega_{U2,z}',...
           '\omega_{L2,x}','\omega_{L2,y}','\omega_{L2,z}') 
    grid;
       
    %plot velocities of 2 centers of mass arm2
    subplot(3,2,5)
    hold on
    plot(time_range,U3vx)
    plot(time_range,U3vy)
    plot(time_range,U3vz)
    plot(time_range,L3vx,'--')
    plot(time_range,L3vy,'--')
    plot(time_range,L3vz,'--')
    title('Centre of masses of arm 3: linear velocities');
    legend('U3v_x','U3v_y','U3v_z','L3v_x','L3v_y','L3v_z')
    subplot(3,2,6)
    hold on
    plot(time_range,omegaU3,'-')
    plot(time_range,omegaL3,'--')
    title('Centre of masses of arm 3: angular velocites');
    legend('\omega_{U3,x}','\omega_{U3,y}','\omega_{U3,z}',...
           '\omega_{L3,x}','\omega_{L3,y}','\omega_{L3,z}')  
    grid;
    %plot 3 velocities end-effector
    figure;
    hold on
    plot(time_range,Pxv)
    plot(time_range,Pyv)
    plot(time_range,Pzv)
    legend(...
        'v_x',...
        'v_y',...
        'v_z'...
    );
    title('end-effector center of Mass')
    grid;
end
