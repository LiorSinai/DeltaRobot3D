function plot_tangents(velocities, accelerations,time_range)
% plot the magnitude of the velocities and the accelerations in the 
% direction of the velocities (tangential accelerations)

% INPUTS 
% velocities = 42xN velocity values
% accelerations = 42xN acceleration values
% time_range = 1xN time values

    % centers of mass arm1
    U1v=velocities(4:6,:);
    U1a=accelerations(4:6,:);
    u=U1v./vecnorm(U1v); % unit vector in the tangential direction
    U1a_tan=sum(u.*U1a); % dot product
    
    omegaU1=velocities(7:9,:);
    alphaU1=accelerations(7:9,:);
    u=omegaU1./vecnorm(omegaU1); % unit vector of axis of rotation
    alphaU1_tan=sum(u.*alphaU1);
    
    L1v=velocities(10:12,:);
    L1a=accelerations(10:12,:);
    u=L1v./vecnorm(L1v); % unit vector in the tangential direction
    L1a_tan=sum(u.*L1a); % dot product
    
    omegaL1=velocities(13:15,:);
    alphaL1=accelerations(13:15,:);
    u=omegaL1./vecnorm(omegaL1); % unit vector of axis of rotation
    alphaL1_tan=sum(u.*alphaL1);

    % second center of mass     
    U2v=velocities(16:18,:);
    U2a=accelerations(16:18,:);
    u=U2v./vecnorm(U2v); % unit vector in the tangential direction
    U2a_tan=sum(u.*U2a); % dot product
    
    omegaU2=velocities(19:21,:);
    alphaU2=accelerations(19:21,:);
    u=omegaU2./vecnorm(omegaU2); % unit vector of axis of rotation
    alphaU2_tan=sum(u.*alphaU2);
    
    L2v=velocities(22:24,:);
    L2a=accelerations(22:24,:);
    u=L2v./vecnorm(L2v); % unit vector in the tangential direction
    L2a_tan=sum(u.*L2a); % dot product

    omegaL2=velocities(25:27,:);
    alphaL2=accelerations(25:27,:);
    u=omegaL2./vecnorm(omegaL2); % unit vector of axis of rotation
    alphaL2_tan=sum(u.*alphaL2);
    
    % third center of mass
    U3v=velocities(28:30,:);
    U3a=accelerations(28:30,:);
    u=U3v./vecnorm(U3v); % unit vector in the tangential direction
    U3a_tan=sum(u.*U3a); % dot product
    
    omegaU3=velocities(31:33,:);
    alphaU3=accelerations(31:33,:);
    u=omegaU3./vecnorm(omegaU3); % unit vector of axis of rotation
    alphaU3_tan=sum(u.*alphaU3);
    
    L3v=velocities(34:36,:);
    L3a=accelerations(34:36,:);
    u=L3v./vecnorm(L3v); % unit vector in the tangential direction
    L3a_tan=sum(u.*L3a); % dot product
    
    omegaL3=velocities(37:39,:);
    alphaL3=accelerations(37:39,:);
    u=omegaL3./vecnorm(omegaL3); % unit vector of axis of rotation
    alphaL3_tan=sum(u.*alphaL3);
        
    %end-effector
    
    Pv = velocities(40:42,:);
    Pa = accelerations(40:42,:);
    u=Pv./vecnorm(Pv); % unit vector in the tangential direction
    Pa_tan=sum(u.*Pa); % dot product
    
    figure('Name','Tangential velocity and accelerations');

    subplot(3,2,1)
    hold on
    plot(time_range,vecnorm(U1v));
    plot(time_range,U1a_tan);
    plot(time_range,vecnorm(omegaU1),'--');
    plot(time_range,alphaU1_tan,'--');
    legend('|U1v|','U1a_{tan}','U1|\omega|','U1\alpha_{tan}')
    title('Upper arm 1');
    grid;

    subplot(3,2,2)
    hold on
    plot(time_range,vecnorm(L1v));
    plot(time_range,L1a_tan);
    plot(time_range,vecnorm(omegaL1),'--');
    plot(time_range,alphaL1_tan,'--');
    legend('|L1v|','L1a_{tan}','L1|\omega|','L1\alpha_{tan}')
    title('Lower arm 1');
    grid;
    
    subplot(3,2,3)
    hold on
    plot(time_range,vecnorm(U2v));
    plot(time_range,U2a_tan);
    plot(time_range,vecnorm(omegaU2),'--');
    plot(time_range,alphaU2_tan,'--');
    legend('|U2v|','U2a_{tan}','U2|\omega|','U2\alpha_{tan}')
    title('Upper arm 2');
    grid;
    
    subplot(3,2,4)
    hold on
    plot(time_range,vecnorm(L2v));
    plot(time_range,L2a_tan);
    plot(time_range,vecnorm(omegaL2),'--');
    plot(time_range,alphaL2_tan,'--');
    legend('|L2v|','L2a_{tan}','L2|\omega|','L2\alpha_{tan}')
    title('Lower arm 2');
    grid;
    
    subplot(3,2,5)
    hold on
    plot(time_range,vecnorm(U3v));
    plot(time_range,U3a_tan);
    plot(time_range,vecnorm(omegaU3),'--');
    plot(time_range,alphaU3_tan,'--');
    legend('|U3v|','U3a_{tan}','U3|\omega|','U3\alpha_{tan}')
    title('Upper arm 3');
    grid;
    
    subplot(3,2,6)
    hold on
    plot(time_range,vecnorm(L3v));
    plot(time_range,L3a_tan);
    plot(time_range,vecnorm(omegaL3),'--');
    plot(time_range,alphaL3_tan,'--');
    legend('|L3v|','L3a_{tan}','L3|\omega|','L3\alpha_{tan}')
    title('Lower arm 3');
    grid;
    
    figure('Name','Tangential velocity and accelerations cont.');
    hold on
    plot(time_range,vecnorm(Pv));
    plot(time_range,Pa_tan);
    legend('|Pv|','Pa_{tan}')
    title('End effector');
    grid;
    

end
