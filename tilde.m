function rtilde =tilde(r)
% implements the tilde operator.
% INPUTS
% r = 3x1 vector
% OUTPUTS
% rtilde=3x3 anti-symmetric matrix. rtilde'=-rtilde

rtilde=[
   0    -r(3) r(2);
   r(3)   0  -r(1);
  -r(2)  r(1)  0];
end