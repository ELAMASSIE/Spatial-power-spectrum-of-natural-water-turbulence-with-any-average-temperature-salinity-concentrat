function [dr] = f_dr(R_rho)
%
dr = R_rho.*1.85 - 0.85;
M1 = find(R_rho>=1); O1 = R_rho(M1);
M2 = find(R_rho<0.5); O2 = R_rho(M2);
dr(M1) = 01 + O1.^0.5 .* (O1 - 1).^0.5;
dr(M2) = O2.*0.15; 
end

