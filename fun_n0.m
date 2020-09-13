function [n0] = fun_n0(T,S,lumbda)
%average refractive-index
a0 = 1.31405;
a1=1.779e-4; a2=-1.05e-6;
a3=1.6e-8; a4=-2.02e-6;
a5 = 15.868;
a6=0.01155; a7=-0.00423;
a8 = -4382; a9 = 1.1455e6;

n0 = a0 + (a1 + a2*T + a3*T^2)*S + a4*T^2 + (a5 + a6*S + a7*T)/lumbda + a8/lumbda^2 + a9/lumbda^3;
end

