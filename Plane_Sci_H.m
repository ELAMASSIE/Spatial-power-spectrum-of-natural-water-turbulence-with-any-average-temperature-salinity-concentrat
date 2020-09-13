function [W] = Plane_Sci_H(av_T, av_S, lumbda, n_L, n_varepsilon,H, n_chi_t)
%This is the function of Scintillation of Plan wave
% where av_T:<T>, av_S:<S>, lumbda:Wave length, n_L:Propagation distance,
% n_varepsilon:\varepsilon, H:temperature-salinity gradient ratio,
% n_chi_t:\Chi_T
%% spectrum
P0 = [21.61; 0.02; 0.61; -18.18; 0.04; 0.55; 174.90; 0.96] ;
f_g = @(c,P,x) (1 + P(1).* c.^P(2) .*x.^P(3) + P(4).*c.^P(5).*x.^P(6)).*exp(-x.^2.*P(7).*c^P(8));
f_g2 =@(Belta,Pr,x,P)f_g(Belta/Pr*0.072^(4/3),P,x);
f_Fay = @(P,A_t, B_s, chi_t, chi_s, chi_ts, varepsilon, Belta, eta, Pr_t, Pr_s, Pr_ts, x) (4*pi)^(-1)*Belta*varepsilon^(-1/3)*eta^(11/3).*x.^(-11/3) .*...
    (A_t^2.*chi_t.*f_g2(Belta,Pr_t,x,P) + B_s^2.*chi_s.*f_g2(Belta,Pr_s,x,P) + 2*A_t*B_s.*chi_ts.*f_g2(Belta,Pr_ts,x,P)); %here x = eta*kappa
%% Scintillation
f_BJ2 = @(P,A_t, B_s, chi_t, chi_s, chi_ts, varepsilon, Belta, eta, Pr_t, Pr_s, Pr_ts, x, L, k)...
    x .* f_Fay(P,A_t, B_s, chi_t, chi_s, chi_ts, varepsilon, Belta, eta, Pr_t, Pr_s, Pr_ts, x)...
    .* ( 1 - sin(L./(k.*eta.^2).*x.^2)./(L./(k.*eta.^2).*x.^2) );
f_J2 = @(P,A_t, B_s, chi_t, chi_s, chi_ts, varepsilon, Belta, eta, Pr_t, Pr_s, Pr_ts, L, k)...
    integral(@(x)f_BJ2(P,A_t, B_s, chi_t, chi_s, chi_ts, varepsilon, Belta, eta, Pr_t, Pr_s, Pr_ts, x, L, k),0,+Inf,'arrayvalued',true);
f_SCI = @(P,A_t, B_s, chi_t, chi_s, chi_ts, varepsilon, Belta, eta, Pr_t, Pr_s, Pr_ts, L, k, n0)...
    8*pi^2.*k.^2.*L./eta.^2./n0.^2.*f_J2(P,A_t, B_s, chi_t, chi_s, chi_ts, varepsilon, Belta, eta, Pr_t, Pr_s, Pr_ts, L, k);
%% parameters
n_Belta = 0.72;

n_k = 2*pi./(lumbda.*1e-9);
[n_A_t,n_B_s] = fun_AB(av_T,av_S,lumbda);
n_Pr_s = Y_Schmidt(av_T,'C',av_S,'ppt');
n_Pr_t = SW_Prandtl(av_T,'C',av_S,'ppt');
n_Pr_ts = 2*n_Pr_t*n_Pr_s/(n_Pr_t+n_Pr_s);

A_Alpha = gsw_alpha(av_S,av_T,0);  
A_Beta = gsw_beta(av_S,av_T,0); 
R_rho = A_Alpha.*abs(H)./A_Beta;
dr = f_dr(R_rho);
n_chi_s = n_chi_t.* dr./ H.^2;
n_chi_ts = (1+dr).*n_chi_t./(2.*H);

n_eta = (SW_Kviscosity(av_T,'C',av_S,'ppt'))^(3/4)/n_varepsilon^(1/4);

n_n0 = fun_n0(av_T,av_S,lumbda);
%% calculate
W = f_SCI(P0,n_A_t, n_B_s, n_chi_t, n_chi_s, n_chi_ts, n_varepsilon, n_Belta, n_eta, n_Pr_t, n_Pr_s, n_Pr_ts, n_L, n_k, n_n0);
end

