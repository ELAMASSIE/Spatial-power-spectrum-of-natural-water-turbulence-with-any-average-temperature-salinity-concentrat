function [A,B] = fun_AB(T,S,lumda)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
a1=1.779e-4; a2=-1.05e-6;
a3=1.6e-8; a4=-2.02e-6;
a6=0.01155; a7=-0.00423;

A = a2.*S + 2.*a3.*T.*S + 2.*a4.*T + a7/lumda;
B = a1 + a2.*T + a3.*T.^2 + a6/lumda;
end

