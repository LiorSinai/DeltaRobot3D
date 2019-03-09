function R=rot3D_Rodrigues(n,theta)
% Lior Sinai, 25 December 2018
% rotates a vector about a unit normal by an angle theta

tol=1e-10; %tolerance
if (~isequal([3 1],size(n)))
    error('expected a unit vector n of size [3 1]'); 
% elseif (~isequal([3 1],size(v)))
%     error('expected a vector of size [3 1]');
% elseif (norm(n)<1-tol ) || (norm(n)>1+tol ) % norm~=1
%     error('expected a unit vector n with magnitude 1'); 
end

%v=v*cos(theta) + cross(n,v)*sin(theta) + (1-cos(theta))*sum(n.*v)*n ; % only works for size(v)=[3 1]

% Equivalent matrix form
nTilde=[0    -n(3) +n(2);
        n(3)   0    -n(1);
       -n(2)  n(1)     0];
% OmegaTilde*v=cross(n,v);
% OmegaTilde^2*v=-v+sum(n.*v)*n=cross(n,cross(n,v)) = v_perp where norm(n)=1
R=eye(3)+ nTilde*sin(theta)+(1-cos(theta))*nTilde^2;

end