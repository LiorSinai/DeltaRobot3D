function plot_acceleration( accelerations,time_range,delta_t )
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
    U1vx = accelerations(4,:);
    U1vy = accelerations(5,:);
    U1vz = accelerations(6,:);
    L1vx = accelerations(10,:);
    L1vy = accelerations(11,:);
    L1vz = accelerations(12,:);

    % second center of mass
    U2vx = accelerations(16,:);
    U2vy = accelerations(17,:);
    U2vz = accelerations(18,:);
    L2vx = accelerations(22,:);
    L2vy = accelerations(23,:);
    L2vz = accelerations(24,:);

    % third center of mass
    U3vx = accelerations(28,:);
    U3vy = accelerations(29,:);
    U3vz = accelerations(30,:);
    L3vx = accelerations(34,:);
    L3vy = accelerations(35,:);
    L3vz = accelerations(36,:);
    
    % actuators angles
    
    omega1 = accelerations(8,:);
    omega2 = accelerations(20,:);
    omega3 = accelerations(32,:);
    
    %end-effector
    
    Pxv = accelerations(40,:);
    Pyv = accelerations(41,:);
    Pzv = accelerations(42,:);
    
    
    % rescale time axis to get proper time
    time = time_range*delta_t;
    
    figure('Name','Accelerations');

    
    %plot 3velocities of 2 centers of mass arm1 
    subplot(3,2,1)
    hold on
    plot(time,U1vx,'-')
    plot(time,U1vy,'-')
    plot(time,U1vz,'--')
    plot(time,L1vx,'--')
    plot(time,L1vy,'--')
    plot(time,L1vz,'--')
    legend('U1a_x','U1a_y','U1a_z','L1a_x','L1a_y','L1a_z')
    title('Centre of masses of arm 1');
    %plot 3velocities of 2 centers of mass arm2
    subplot(3,2,2)
    hold on
    plot(time,U2vx)
    plot(time,U2vy)
    plot(time,U2vz)
    plot(time,L2vx,'--')
    plot(time,L2vy,'--')
    plot(time,L2vz,'--')
    title('Centre of masses of arm 2');
    legend('U2a_x','U2a_y','U2a_z','L2a_x','L2a_y','L2a_z')
    
    %plot 3velocities of 2 centers of mass arm3 
    subplot(3,2,3)
    hold on
    plot(time,U3vx)
    plot(time,U3vy)
    plot(time,U3vz)
    plot(time,L3vx,'--')
    plot(time,L3vy,'--')
    plot(time,L3vz,'--')
    title('Centre of masses of arm 2');
    legend('U3a_x','U3a_y','U3a_z','L3a_x','L3a_y','L3a_z')
    
    %plot velocities of actuators angles
    subplot(3,2,4)
    hold on
    plot(time,omega1)
    plot(time,omega2)
    plot(time,omega3)
    legend(...
        '\alpha arm 1',...
        '\alpha arm 2',...
        '\alpha arm 3'...
    );
    title('Actuator angular accelerations')
    
    %plot 3 velocities end-effector
    subplot(3,2,5)
    hold on
    plot(time,Pxv)
    plot(time,Pyv)
    plot(time,Pzv)
    legend(...
        'a_x',...
        'a_y',...
        'a_z'...
    );
    title('end-effector center of Mass')
end
