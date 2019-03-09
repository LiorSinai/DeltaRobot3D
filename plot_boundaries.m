%% Plot boundaries
[X,Y,Z]=sphere(16);
X=(L_u+L_l)*X;
Y=(L_u+L_l)*Y;
Z=(L_u+L_l)*Z;

x10=0.5*(L_b-L_e)*cos(thetaU1_z_0);
y10=0.5*(L_b-L_e)*sin(thetaU1_z_0);
hSphere1=surf(X+x10,Y+y10,Z);
hold on;
set(hSphere1, 'FaceAlpha', 0.1,'FaceColor','k'); %red

x20=0.5*(L_b-L_e)*cos(thetaU2_z_0);
y20=0.5*(L_b-L_e)*sin(thetaU2_z_0);
hSphere2=surf(X+x20,Y+y20,Z);
set(hSphere2, 'FaceAlpha', 0.1,'FaceColor','m');

x30=0.5*(L_b-L_e)*cos(thetaU3_z_0);
y30=0.5*(L_b-L_e)*sin(thetaU3_z_0);
hSphere3=surf(X+x30,Y+y30,Z);
set(hSphere3, 'FaceAlpha', 0.1,'FaceColor','r'); 
% lighting gouraud
% %shading interp

set(hSphere1,'EdgeColor','None');
set(hSphere2,'EdgeColor','None');
set(hSphere3,'EdgeColor','None');

% plot sphere centres
plot(x10,y10,'ok',x20,y20,'om',x30,y30,'or')
