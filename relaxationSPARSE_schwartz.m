function [phi, Dphi] = relaxationSPARSE_schwartz(a,b,t)

% For two vectors a, b of the same length and a positive scalar t, this 
% function computes
    % phi(a,b) = (a-t)*(b-t)               if    a+b >= 2*t
    %            -0.5*((a-t)^2 - (b-t)^2)  else
% and the gradient Dphi = [Dphi_a Dphi_b]

phi =   (a-t).*(b-t).*(a+b >= 2*t) ...
      - 0.5*((a-t).^2 + (b-t).^2).*(a+b < 2*t) ;
      
if nargout > 1
    % (n_ab x 2) gradient of the relaxation function
    Dphi =   [b-t a-t].*repmat((a+b >= 2*t),1,2) ...
           - [a-t b-t].*repmat((a+b <  2*t),1,2);
end