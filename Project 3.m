%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Multivariable Control
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design a reduced order model and state-feedback to minimize the oscillations 
% of the original system while maintaining similar performance response 
% (setting time and output magnitude), see Figure 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

% Using 1st input 
h11 = tf([8 27.5 30.25 10.5],[1 4 6.25 5.25 2.25]);
h21 = tf([5 22.5 35.5 19.5],[1 4 6.25 5.25 2.25]);
h31 = tf([3 13 22.25 14.25],[1 4 6.25 5.25 2.25]);
% Using 2nd input 
h12 = tf([5 10.5 1.5 -6.25],[1 4 6.25 5.25 2.25]);
h22 = tf([10 32 34 11],[1 4 6.25 5.25 2.25]);
h32 = tf([6 19.5 23.75 8.5],[1 4 6.25 5.25 2.25]);
% Combined transfer function of MIMO system
H = [h11 h12;h21 h22;h31 h32];

% State space representation of MIMO system
sys = ss(H);
sys = minreal(sys);

[A,B,C,D] = ssdata(sys);
[Modal,Diagonal] = eig(A);

%% ------- Residue Approach Model Reduction (Second Order System) ------- %%
% Make Diagonal Matrix in Descending Order
M_hat = [Modal(:,1),Modal(:,2),Modal(:,3),Modal(:,4)];
 
Ad = inv(M_hat)*A*M_hat;
Bd = inv(M_hat)*B;
Cd = C*M_hat;

A11 = Ad(1:2,1:2);
A12 = Ad(1:2,3:4);
A21 = Ad(3:4,1:2);
A22 = Ad(3:4,3:4);

B1  = Bd(1:2,:);
B2  = Bd(3:4,:);

C1  = Cd(:,1:2);
C2  = Cd(:,3:4);

% Calculating Reduced Order Matrices
Ar  = A11 - A12*inv(A22)*A21;
Br  = B1 - A12*inv(A22)*B2;
Cr  = C1 - C2*inv(A22)*A21;
Dr  = D - C2*inv(A22)*B2;

% State-Space Representation of Reduced Order System
sys_r = ss(Ar,Br,Cr,Dr);

% Taking real part of eigenvalues of Reduced Order System
P = [-0.5,-0.5];
K = place(Ar,Br,P); 

A_cl = Ar - Br*K;
B_cl = Br;
C_cl = Cr;
% Adjusting the steady state error to maintain performance of original system
D_cl = C*inv(-A)*B+D-C_cl*inv(-A_cl)*B_cl; 
sys_cl = ss(A_cl,B_cl,C_cl,D_cl);

figure;
step(sys);
hold on
step(sys_r);
hold on
step(sys_cl);
title('Step Response of Original, Reduced and State-Feedback Second Order');

figure;
impulse(sys);
hold on
impulse(sys_r);
hold on
impulse(sys_cl);
title('Step Response of Original, Reduced and State-Feedback Second Order');

