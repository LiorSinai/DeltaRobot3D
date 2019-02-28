%% Inverse Kinematics
%Lior Sinai, 2019-02-24

Nvalid=size(thetaA_t,2); % number of valid points
time_valid=time_range(1:Nvalid); % if a singularity was hit, the valid_time
                                 % will be shorter than the time_range 
%                                 
testIndex=100;
%qa0=[-0.29;-0.73;-0.29]; %random initial guess

indP=40:42;
%% Test position
qa_t=zeros(size(thetaA_t));
qaVel_t=qa_t;
qaAcc_t=qa_t;
for testIndex=1:Nvalid
    if testIndex==1
        qa0=thetaA_0;
    else
        qa0=thetaA_t(:,testIndex-1);
    end
    
    [qa, qaVel, qaAcc]=calculate_IK(L,qa0,Q(indP,testIndex),...
        Qvel(indP,testIndex),Qacc(indP,testIndex),...
        tolerance);
    qa_t(:,testIndex)=qa;
    qaVel_t(:,testIndex)=qaVel;
    qaAcc_t(:,testIndex)=qaAcc;
end

%% Plot results
figure;
plot(time_valid,qa_t-thetaA_t,'x-')
title('Error in \theta_A')

figure;
%plot(time_valid,qaVel_t-alphaA*time_range,'x-') % constant acceleratio
%plot(time_valid,qaVel_t-omega,'x-') % constant velocity
% omegaDrive=[(+30*pi/180)*cos(time_valid);
%            (-20*pi/180)*cos(time_valid);
%            (+45*pi/180)*cos(time_valid)];
thetaFunc=[driveFunc.thetaA1(t_sym); driveFunc.thetaA2(t_sym); driveFunc.thetaA3(t_sym)];
omegaFunc=diff(thetaFunc,t_sym);
omegaDrive=eval(subs(omegaFunc,t_sym,time_valid));
plot(time_valid,qaVel_t-omegaDrive)
title('Error in \omega_A')

figure;
%plot(time_valid,qaAcc_t-alphaA,'x-') % cosntant acceleration
%plot(time_valid,qaAcc_t-0,'x-') % constant velocity
% alphaDrive=[-(+30*pi/180)*sin(time_valid);
%            -(-20*pi/180)*sin(time_valid);
%            -(+45*pi/180)*sin(time_valid)];
alphaFunc=diff(omegaFunc,t_sym);
alphaDrive=eval(subs(alphaFunc,t_sym,time_valid));
plot(time_valid,qaAcc_t-alphaDrive)
title('Error in \alpha_A')
