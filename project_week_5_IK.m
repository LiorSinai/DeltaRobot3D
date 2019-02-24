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
        omega,Qvel(indP,testIndex),Qacc(indP,testIndex),...
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
plot(time_valid,qaVel_t-omega,'x-')
title('Error in \omega_A')

figure;
plot(time_valid,qaAcc_t,'x-')
title('Error in \alpha_A')
