function gamma_IK = Gamma_IKNumeric(L_b,L_e,L_u,P_x,P_y,P_z,Pacc_x,Pacc_y,Pacc_z,Pvel_x,Pvel_y,Pvel_z,omegaA1,omegaA2,omegaA3,thetaA1,thetaA2,thetaA3)
%GAMMA_IKNUMERIC
%    GAMMA_IK = GAMMA_IKNUMERIC(L_B,L_E,L_U,P_X,P_Y,P_Z,PACC_X,PACC_Y,PACC_Z,PVEL_X,PVEL_Y,PVEL_Z,OMEGAA1,OMEGAA2,OMEGAA3,THETAA1,THETAA2,THETAA3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    24-Feb-2019 12:44:18

t2 = cos(thetaA1);
t3 = L_u.^2;
t4 = sin(thetaA1);
t5 = sqrt(3.0);
t6 = L_b./2.0;
t7 = L_e./2.0;
t8 = L_u.*t4;
t9 = -t6+t7+t8;
t10 = P_z.*2.0;
t11 = P_x.*2.0;
t12 = Pvel_x.^2;
t13 = Pvel_y.^2;
t14 = Pvel_z.^2;
t15 = cos(thetaA2);
t16 = sin(thetaA2);
t17 = L_b./4.0;
t18 = L_e./4.0;
t19 = P_y.*2.0;
t20 = L_u.*t16;
t21 = -t6+t7+t20;
t22 = cos(thetaA3);
t23 = sin(thetaA3);
gamma_IK = [t12.*-2.0-t13.*2.0-t14.*2.0-omegaA1.^2.*(t2.^2.*t3.*2.0+t3.*t4.^2.*2.0+L_u.*t4.*(P_x+t17-t18-(L_u.*t4)./2.0)-L_u.*t2.*(P_z+L_u.*t2).*2.0+L_u.*t4.*t5.*(P_y-(t5.*t9)./2.0))-Pacc_z.*(t10+L_u.*t2.*2.0)-Pacc_y.*(t19-t5.*t9)-Pacc_x.*(t6-t7+t11-L_u.*t4)+L_u.*Pvel_x.*omegaA1.*t2.*2.0+L_u.*Pvel_z.*omegaA1.*t4.*4.0+L_u.*Pvel_y.*omegaA1.*t2.*t5.*2.0;t12.*-2.0-t13.*2.0-t14.*2.0-omegaA2.^2.*(t3.*t15.^2.*2.0+t3.*t16.^2.*2.0+L_u.*t16.*(P_x+t17-t18-(L_u.*t16)./2.0)-L_u.*t15.*(P_z+L_u.*t15).*2.0-L_u.*t5.*t16.*(P_y+(t5.*t21)./2.0))-Pacc_z.*(t10+L_u.*t15.*2.0)-Pacc_y.*(t19+t5.*t21)-Pacc_x.*(t6-t7+t11-L_u.*t16)+L_u.*Pvel_x.*omegaA2.*t15.*2.0+L_u.*Pvel_z.*omegaA2.*t16.*4.0-L_u.*Pvel_y.*omegaA2.*t5.*t15.*2.0;t12.*-2.0-t13.*2.0-t14.*2.0-P_y.*Pacc_y.*2.0-omegaA3.^2.*(t3.*t22.^2.*2.0+t3.*t23.^2.*2.0-L_u.*t23.*(P_x-t6+t7+L_u.*t23).*2.0-L_u.*t22.*(P_z+L_u.*t22).*2.0)-Pacc_x.*(-L_b+L_e+t11+L_u.*t23.*2.0)-Pacc_z.*(t10+L_u.*t22.*2.0)-L_u.*Pvel_x.*omegaA3.*t22.*4.0+L_u.*Pvel_z.*omegaA3.*t23.*4.0];
