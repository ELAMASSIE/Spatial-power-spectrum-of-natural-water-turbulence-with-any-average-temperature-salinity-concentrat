clc; clear;
close all;
%% parameters
n_Yip   = 1E-2; 
n_Belta = 0.72; 
n_omiga = -0.25;
n_Hay_t = 1e-5;
n_H = -20;
a_L = (0.1:0.1:220);

%% for av_T = 15; av_S = (0,10,20,30,40);lumbda = 532;
lumbda = 532;
av_T = 15; a_av_S = (0:10:40);
Data1 = double(zeros(length(a_L),length(a_av_S)));
for j = 1:length(a_av_S)
    Data1(:,j) = Plane_Sci_H(av_T, a_av_S(j), lumbda, a_L, n_Yip,n_H, n_Hay_t);
end
figure(1);
hold on; 
plot(a_L',Data1,'-','linewidth',2);plot([0,max(a_L)],[1,1],'--k');
legend('0 ppt','10 ppt','20 ppt','30 ppt','40 ppt');
axis([0,150,0,1.3]);
xlabel('L'),ylabel('\sigma _{I,{\rm{pl}}}^2');hold off;
SD1 =[a_L',Data1];
% save('f5_a.mat','SD1');
%% for av_T = (0,10,20,30); av_S = 34.9;lumbda = 532;
lumbda = 532;
a_av_T = (0:10:30); av_S = 34.9; 
Data2 = double(zeros(length(a_L),length(a_av_T)));
for j = 1:length(a_av_T)
    Data2(:,j) = Plane_Sci_H(a_av_T(j), av_S, lumbda, a_L, n_Yip,n_H, n_Hay_t);
end
figure(2);
hold on; 
plot(a_L,Data2,'-','linewidth',2);plot([0,max(a_L)],[1,1],'--k');
legend('0^\circ C','10^\circ C','20^\circ C','30^\circ C');
axis([0,220,0,1.3]);
xlabel('L'),ylabel('\sigma _{I,{\rm{pl}}}^2');hold off;
SD2 =[a_L',Data2];
% save('f5_b.mat','SD2');

%% for av_T = 15; av_S = 34.9;a_lumbda = (450:50:600);
a_lumbda = (450:50:600);
av_T = 15; av_S = 34.9; 
Data3 = double(zeros(length(a_L),length(a_lumbda)));
for j = 1:length(a_lumbda)
    Data3(:,j) = Plane_Sci_H(av_T, av_S, a_lumbda(j), a_L, n_Yip,n_H, n_Hay_t);
end
figure(3);
hold on; 
plot(a_L,Data3,'-','linewidth',2);plot([0,max(a_L)],[1,1],'--k');
legend('450nm','500nm','550nm','600nm');
axis([0,150,0,1.3]);
xlabel('L'),ylabel('\sigma _{I,{\rm{pl}}}^2');hold off;
SD3 =[a_L',Data3];
% save('f5_c.mat','SD3');
