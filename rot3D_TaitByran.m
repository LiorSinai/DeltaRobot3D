% Lior Sinai 
% 24 December 2018

function R=rot3D_TaitByran(phi, theta, psi)
% rotate a 3D vector using a Tait-Byran Euler angle sequence
% Convention used: ZYX - psi, theta, phi
%                        roll, pitch, yaw
%                        bank, attitude, heading
% See Representing Attitude: Euler Angles, Unit Quaternions, and Rotation
% Vectors, James Diebel (2006)
% Rotation matrix is for a co-ordinate transformation. Vectors use R'

% INPUTS
% phi, theta, psi = Euler angles

% OUTPUTS
% R = 3x3 rotation matrix

% For symoblic checks
%syms phi theta psi

Rx=[1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
Ry=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz=[cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];

R=Rz*Ry*Rx;
end

% test
% vector=[0 1 1; 0 0 1; 0 0 0];
% plot3(vector(1,:),vector(2,:),vector(3,:))
% vectorRot=rot3D_TaitByran(vector,1,0,0);
% plot3(vectorRot(1,:),vectorRot(2,:),vectorRot(3,:))
