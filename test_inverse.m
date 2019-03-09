%% Test inverse kinematics
qa_t=zeros(size(thetaA_t));
qaVel_t=qa_t;
qaAcc_t=qa_t;
indP=40:42;
N=length(time_range);
for testIndex=1:N
    qa0=thetaA_t(:,testIndex);
    [qa, qaVel, qaAcc]=calculate_IK(L,qa0,Q(indP,testIndex),...
        Qvel(indP,testIndex),Qacc(indP,testIndex),...
        tolerance);
    qa_t(:,testIndex)=qa;
    qaVel_t(:,testIndex)=qaVel;
    qaAcc_t(:,testIndex)=qaAcc;
end

figure('Name','Test: inverse calculations');
plot(1:N,qa_t-thetaA_t,'x-')
title('Error in \theta_A')

% figure;
% plot(1:N,qaVel_t-omega,'x-')
% title('Error in \omega_A')

% figure;
% plot(1:N,qaAcc_t,'x-')
% title('Error in \alpha_A')
