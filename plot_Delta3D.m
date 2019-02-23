function plot_Delta3D( Q,L,R_A1,R_A2,R_A3,R_sym,thetaA,R_L1,R_L2,R_L3,thetaA_t,time_range,delay_between_plots )
    
%%%             z
%%%             |
%%% Start at 1__M__2_ _ _x          y
%%%          /     \U           3\ |
%%%          \     /L             \|_ _2 x
%%%           \_P_/               /
%%%                             1/

    % Centers of Mass
    M_x = Q(1,:);
    M_y = Q(2,:);
    M_z = Q(3,:);
    %arm1
    U1_x = Q(4,:);
    U1_y = Q(5,:);
    U1_z = Q(6,:);
    L1_x = Q(10,:);
    L1_y = Q(11,:);
    L1_z = Q(12,:);
    %arm2
    U2_x = Q(16,:);
    U2_y = Q(17,:);
    U2_z = Q(18,:);
    L2_x = Q(22,:);
    L2_y = Q(23,:);
    L2_z = Q(24,:);
    %arm3
    U3_x = Q(28,:);
    U3_y = Q(29,:);
    U3_z = Q(30,:);
    L3_x = Q(34,:);
    L3_y = Q(35,:);
    L3_z = Q(36,:);
    %end effector
    P_x = Q(40,:);
    P_y = Q(41,:);
    P_z = Q(42,:);

    R_A1_calc = R_A1*[L(4)/2; 0 ;0];
    R_U1_calc = zeros(3,size(Q,2));
    R_L1_calc = zeros(3,size(Q,2));
    R_A2_calc = R_A2*[L(4)/2; 0 ;0];
    R_U2_calc = zeros(3,size(Q,2));
    R_L2_calc = zeros(3,size(Q,2));
    R_A3_calc = R_A3*[L(4)/2; 0 ;0];
    R_U3_calc = zeros(3,size(Q,2));
    R_L3_calc = zeros(3,size(Q,2));
    

    for timeStep = 1:length(time_range)
        R_U1 =R_A1*rot3D_Rodrigues([0;1;0],thetaA_t(1,timeStep));         
        R_U2 =R_A2*rot3D_Rodrigues([0;1;0],thetaA_t(2,timeStep)); 
        R_U3 =R_A3*rot3D_Rodrigues([0;1;0],thetaA_t(3,timeStep)); 
        
        R_U1_calc(:,timeStep) = R_U1*[0;0;-L(1)/2];
        R_U2_calc(:,timeStep) = R_U2*[0;0;-L(1)/2];
        R_U3_calc(:,timeStep) = R_U3*[0;0;-L(1)/2];
        R_L1_calc(:,timeStep) = R_L1(:,:,timeStep)*[0;0;-L(2)/2];
        R_L2_calc(:,timeStep) = R_L2(:,:,timeStep)*[0;0;-L(2)/2];
        R_L3_calc(:,timeStep) = R_L3(:,:,timeStep)*[0;0;-L(2)/2];

    end
    % Coordinates Ai,Bi,Ci   
%%%             z
%%%             |
%%%       A_ 1__M__2_ _ _x          y
%%%       B_ /     \U           3\ |
%%%          \     /L             \|_ _2 x
%%%         C_\_P_/               /
%%%                             1/    
    
    %arm1
    A1_x = M_x+R_A1_calc(1,1);
    A1_y = M_y+R_A1_calc(2,1);
    A1_z = M_z+R_A1_calc(3,1);
    B1_x = U1_x+R_U1_calc(1,:);
    B1_y = U1_y+R_U1_calc(2,:);
    B1_z = U1_z+R_U1_calc(3,:);
    C1_x = L1_x+R_L1_calc(1,:);
    C1_y = L1_y+R_L1_calc(2,:);
    C1_z = L1_z+R_L1_calc(3,:);
    %arm2
    A2_x = M_x+R_A2_calc(1,1);
    A2_y = M_y+R_A2_calc(2,1);
    A2_z = M_z+R_A2_calc(3,1);
    B2_x = U2_x+R_U2_calc(1,:);
    B2_y = U2_y+R_U2_calc(2,:);
    B2_z = U2_z+R_U2_calc(3,:);
    C2_x = L2_x+R_L2_calc(1,:);
    C2_y = L2_y+R_L2_calc(2,:);
    C2_z = L2_z+R_L2_calc(3,:);
    %arm3
    A3_x = M_x+R_A3_calc(1,1);
    A3_y = M_y+R_A3_calc(2,1);
    A3_z = M_z+R_A3_calc(3,1);
    B3_x = U3_x+R_U3_calc(1,:);
    B3_y = U3_y+R_U3_calc(2,:);
    B3_z = U3_z+R_U3_calc(3,:);
    C3_x = L3_x+R_L3_calc(1,:);
    C3_y = L3_y+R_L3_calc(2,:);
    C3_z = L3_z+R_L3_calc(3,:);
       
    figure('Name','Position');
    for timeStep = 1:length(time_range)
        if timeStep==1
        % Make a brand new plot
        h.M=plot3(...
            ... % centre of masses
            M_x(timeStep),M_y(timeStep),M_z(timeStep),'mx',...
            ... % base
            [A1_x(timeStep), A2_x(timeStep), A3_x(timeStep),A1_x(timeStep)],...
            [A1_y(timeStep), A2_y(timeStep), A3_y(timeStep),A1_y(timeStep)],...
            [A1_z(timeStep), A2_z(timeStep), A3_z(timeStep),A1_z(timeStep)],'k'...
        );
        hold on;
        h.P=plot3(...
            ... % centre of masses
            P_x(timeStep),P_y(timeStep),P_z(timeStep),'kx',...
            ... % end effector base
            [C1_x(timeStep), C2_x(timeStep), C3_x(timeStep),C1_x(timeStep)],...
            [C1_y(timeStep), C2_y(timeStep), C3_y(timeStep),C1_y(timeStep)],...
            [C1_z(timeStep), C2_z(timeStep), C3_z(timeStep),C1_z(timeStep)],'k'...
        );
        h.arm1=plot3(... 
            ...% arm 1: centre of masses    
            U1_x(timeStep),U1_y(timeStep),U1_z(timeStep),'kx',...
            L1_x(timeStep),L1_y(timeStep),L1_z(timeStep),'k*',...
            ... % arm 1
            [0,A1_x(timeStep),B1_x(timeStep),C1_x(timeStep),P_x(timeStep)],...
            [0,A1_y(timeStep),B1_y(timeStep),C1_y(timeStep),P_y(timeStep)],...
            [0,A1_z(timeStep),B1_z(timeStep),C1_z(timeStep),P_z(timeStep)],'blue'...
           );
        h.arm2=plot3(...
            ...% arm 2: centre of masses
            U2_x(timeStep),U2_y(timeStep),U2_z(timeStep),'mx',...
            L2_x(timeStep),L2_y(timeStep),L2_z(timeStep),'m*',...
            ... % arm 2
            [0,A2_x(timeStep),B2_x(timeStep),C2_x(timeStep),P_x(timeStep)],...
            [0,A2_y(timeStep),B2_y(timeStep),C2_y(timeStep),P_y(timeStep)],...
            [0,A2_z(timeStep),B2_z(timeStep),C2_z(timeStep),P_z(timeStep)],'blue'...
            );
        h.arm3=plot3(...
            ... arm 3: centre of masses
            U3_x(timeStep),U3_y(timeStep),U3_z(timeStep),'rx',...
            L3_x(timeStep),L3_y(timeStep),L3_z(timeStep),'r*',...
            ... % arm 3
            [0,A3_x(timeStep),B3_x(timeStep),C3_x(timeStep),P_x(timeStep)],...
            [0,A3_y(timeStep),B3_y(timeStep),C3_y(timeStep),P_y(timeStep)],...
            [0,A3_z(timeStep),B3_z(timeStep),C3_z(timeStep),P_z(timeStep)],'blue'...
          );

        xmax=abs(L(1))+abs(L(2));
        ymax=xmax;
        zmax=xmax;
        axis([-xmax,xmax,-ymax,ymax,-zmax,zmax]);
        % plot current stime step on plot
        message=text(xmax,ymax,zmax*0.8,sprintf('%.1f s, step %d',time_range(timeStep),timeStep));
        xlabel('x');
        ylabel('y');
        zlabel('z');
        grid;
        
       
        temptext=text(0,0,zmax*0.4,'Press any key to continue');
        pause;
        delete(temptext);
        else
        % Keep labels, title, orientation, axes, etc
        % Update only values
        % This way you can rotate the plot while it is plotting 
            h.M(1).XData=M_x(timeStep);
            h.M(1).YData=M_y(timeStep);
            h.M(1).ZData=M_z(timeStep);
        % end effector         
            h.P(1).XData=P_x(timeStep);
            h.P(1).YData=P_y(timeStep);
            h.P(1).ZData=P_z(timeStep);
            
            h.P(2).XData=[C1_x(timeStep), C2_x(timeStep), C3_x(timeStep),C1_x(timeStep)];
            h.P(2).YData=[C1_y(timeStep), C2_y(timeStep), C3_y(timeStep),C1_y(timeStep)];
            h.P(2).ZData=[C1_z(timeStep), C2_z(timeStep), C3_z(timeStep),C1_z(timeStep)];
         % Arm 1  
            h.arm1(1).XData=U1_x(timeStep);
            h.arm1(1).YData=U1_y(timeStep);
            h.arm1(1).ZData=U1_z(timeStep);
         
            h.arm1(2).XData=L1_x(timeStep);
            h.arm1(2).YData=L1_y(timeStep);
            h.arm1(2).ZData=L1_z(timeStep); 

            h.arm1(3).XData=[0,A1_x(timeStep),B1_x(timeStep),C1_x(timeStep),P_x(timeStep)];
            h.arm1(3).YData=[0,A1_y(timeStep),B1_y(timeStep),C1_y(timeStep),P_y(timeStep)];
            h.arm1(3).ZData=[0,A1_z(timeStep),B1_z(timeStep),C1_z(timeStep),P_z(timeStep)];
      % Arm 2
            h.arm2(1).XData=U2_x(timeStep);
            h.arm2(1).YData=U2_y(timeStep);
            h.arm2(1).ZData=U2_z(timeStep);
            
            h.arm2(2).XData=L2_x(timeStep);
            h.arm2(2).YData=L2_y(timeStep);
            h.arm2(2).ZData=L2_z(timeStep);
            
            h.arm2(3).XData=[0,A2_x(timeStep),B2_x(timeStep),C2_x(timeStep),P_x(timeStep)];
            h.arm2(3).YData=[0,A2_y(timeStep),B2_y(timeStep),C2_y(timeStep),P_y(timeStep)];
            h.arm2(3).ZData= [0,A2_z(timeStep),B2_z(timeStep),C2_z(timeStep),P_z(timeStep)];
     % Arm 3: centre of masses        
            h.arm3(1).XData=U3_x(timeStep);
            h.arm3(1).YData=U3_y(timeStep);
            h.arm3(1).ZData=U3_z(timeStep);
            
            h.arm3(2).XData=L3_x(timeStep);
            h.arm3(2).YData=L3_y(timeStep);
            h.arm3(2).ZData=L3_z(timeStep);
      
            h.arm3(3).XData=[0,A3_x(timeStep),B3_x(timeStep),C3_x(timeStep),P_x(timeStep)];
            h.arm3(3).YData=[0,A3_y(timeStep),B3_y(timeStep),C3_y(timeStep),P_y(timeStep)];
            h.arm3(3).ZData=[0,A3_z(timeStep),B3_z(timeStep),C3_z(timeStep),P_z(timeStep)];
           
            message.String=sprintf('%.1f s, step %d',time_range(timeStep),timeStep);     
        end
        drawnow % update the plot so it shows the new values
        pause(delay_between_plots)
    end
end
