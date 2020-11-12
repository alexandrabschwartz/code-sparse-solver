function [phi, Dphi] = relaxationSPARSE(a, b, t, relaxation)

% For two vectors a and b of the same length and a positive scalar t this 
% function evaluates the specified relaxation function.
    
a = a(:);
b = b(:);

switch lower(relaxation)
    case 'scholtes'
        [phi, Dphi] = relaxationSPARSE_scholtes(a,b,t);
    case 'steffensen'
        [phi, Dphi] = relaxationSPARSE_steffensen(a,b,t);
    case 'kadrani'
        [phi, Dphi] = relaxationSPARSE_kadrani(a,b,t);
    case 'schwartz'
        [phi, Dphi] = relaxationSPARSE_schwartz(a,b,t);
end