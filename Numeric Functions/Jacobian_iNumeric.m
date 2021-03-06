function Jacobian_i = Jacobian_iNumeric(L_u,thetaA1,thetaA2,thetaA3)
%JACOBIAN_INUMERIC
%    JACOBIAN_I = JACOBIAN_INUMERIC(L_U,THETAA1,THETAA2,THETAA3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    17-Feb-2019 19:07:39

t2 = cos(thetaA1);
t3 = sin(thetaA1);
t4 = sqrt(3.0);
t5 = (L_u.*t3)./4.0;
t6 = (L_u.*t2)./2.0;
t7 = (L_u.*t3.*t4)./4.0;
t8 = cos(thetaA2);
t9 = sin(thetaA2);
t10 = (L_u.*t4.*t9)./4.0;
t11 = (L_u.*t9)./4.0;
t12 = (L_u.*t8)./2.0;
t13 = cos(thetaA3);
t14 = sin(thetaA3);
t15 = (L_u.*t13)./2.0;
t16 = (L_u.*t14)./2.0;
Jacobian_i = reshape([0.0,t6,t7,0.0,t6,t7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,L_u.*t2.*(-1.0./2.0),0.0,-t5,-t6,0.0,-t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,L_u.*t3.*t4.*(-1.0./4.0),t5,0.0,-t7,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,-t10,0.0,t12,-t10,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,L_u.*t8.*(-1.0./2.0),0.0,-t11,-t12,0.0,-t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t10,t11,0.0,t10,t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t15,0.0,0.0,t15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,L_u.*t13.*(-1.0./2.0),0.0,t16,-t15,0.0,t16,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,L_u.*t14.*(-1.0./2.0),0.0,0.0,-t16,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[33,9]);
