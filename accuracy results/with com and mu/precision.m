clc; clear all; close all;

q_10 = importdata('qData_10.txt');
q_30 = importdata('qData_30.txt');
q_60 = importdata('qData_60.txt');
q_100 = importdata('qData_100.txt');

error_10 = sum(abs(q_100-q_10));
error_30 = sum(abs(q_100-q_30));
error_60 = sum(abs(q_100-q_60));

display(sum(error_10),'Cummulative unaccuracy with 10 points = ');
display(sum(error_30),'Cummulative unaccuracy with 30 points = ');
display(sum(error_60),'Cummulative unaccuracy with 60 points = ');
% 
% s = 1:100;
% 
% figure();
% plot(s,error_10,'*b'); hold on;
% plot(s,error_30,'sr'); hold on;
% plot(s,error_60,'xg'); hold off;
% 
% %--------------------------------
% 
% q_10 = importdata('dqData_10.txt');
% q_30 = importdata('dqData_30.txt');
% q_60 = importdata('dqData_60.txt');
% q_100 = importdata('dqData_100.txt');
% 
% error_10 = sum(abs(q_100-q_10));
% error_30 = sum(abs(q_100-q_30));
% error_60 = sum(abs(q_100-q_60));
% 
% s = 1:100;
% 
% figure();
% plot(s,error_10,'*b'); hold on;
% plot(s,error_30,'sr'); hold on;
% plot(s,error_60,'xg'); hold off;
% 
% %--------------------------------
% 
% q_10 = importdata('ddqData_10.txt');
% q_30 = importdata('ddqData_30.txt');
% q_60 = importdata('ddqData_60.txt');
% q_100 = importdata('ddqData_100.txt');
% 
% error_10 = sum(abs(q_100-q_10));
% error_30 = sum(abs(q_100-q_30));
% error_60 = sum(abs(q_100-q_60));
% 
% s = 1:100;
% 
% figure();
% plot(s,error_10,'*b'); hold on;
% plot(s,error_30,'sr'); hold on;
% plot(s,error_60,'xg'); hold off;
% 
% %--------------------------------