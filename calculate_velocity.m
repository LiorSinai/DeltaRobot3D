function velocities = calculate_velocity( Velocity,Q_sym,R_sym,thetaA_sym,t_sym,Q,R_L1,R_L2,R_L3,thetaA,time_range)
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