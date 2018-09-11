%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Multivariable Control
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Reduce the order of the system to n = 2 and stabilize the output by mirroring the
% unstable eigenvalues about the jw-axis onto the stable region on the s-plane.
% (b) Implement the reduced order model as a closed-loop observer with observer eigenvalues
% 5x the magnitude of the closed-loop system
% (c) Plot the response of the system designed system for (i) impulse inputs, (ii) step inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
 
A = [1.38   -0.2077 6.715 -5.676;
    -0.5814 -4.29   0      0.6750;
     1.067   4.273 -6.654  5.893;
     0.048   4.273  1.343 -2.104];
B = [0      0;
     5.679  0;
     1.136 -3.146;
     1.136  0];
C = [1 0 1  0;
     0 1 0 -1];

% State-space representation of original system
sys = ss(A,B,C,0);
[A,B,C,D] = ssdata(sys)
[Modal,Diagonal] = eig(A)

% -- Residue Approach Model Reduction (Second Order System) -- %
% Make Diagonal Matrix in Descending Order
M_hat = [Modal(:,1),Modal(:,2),Modal(:,4),Modal(:,3)]

Ad = inv(M_hat)*A*M_hat
Bd = inv(M_hat)*B
Cd = C*M_hat

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
 
% State space model for the reduced order system
sys_r = ss(Ar,Br,Cr,Dr)

% % Using direct function reduce()
% dir1 = reduce(sys,2)
% 
% figure;
% step(sys_r-dir1)
 
% State feedback 
P =[-1.9910 -0.0635] % Mirroring the unstable eigenvalues
K = place(Ar,Br,P)

% System with state feedback (K)
Acon = Ar - Br*K;
Bcon = Br;
Ccon = Cr;
Dcon = -C*inv(A)*B + Cr*inv(Ar)*Br;

% State-space representation of controller
sysCon = ss(Acon,Bcon,Ccon,Dcon)

% Observer gain
Po = 5*P;
L = (place(Ar',Cr',Po))';

% Observer system
Aobs = Ar - L*Cr;
Bobs = Br;
Cobs = Cr;
Dobs = -C*inv(A)*B + Cr*inv(Ar)*Br;
 
% State-space representation of observer
sysObs = ss(Aobs,Bobs,Cobs,Dobs)
 
% Close-loop system with observer
Acl = [Ar   -Br*K;
       L*Cr  Ar-L*Cr-Br*K];
Bcl = [Br;
       Br];
Ccl = [Cr zeros(size(Cr))];

% State-space representation of closed-loop system
sys_cl = ss(Acl,Bcl,Ccl,0)
 
% Impulse response
figure;
title('Impulse response')
subplot(2,2,1)
impulse(sys,'r');
title('Original system')
subplot(2,2,2)
impulse(sysCon,'g');
title('System with controller')
subplot(2,2,3)
impulse(sysObs,'b')
title('System with observer')
subplot(2,2,4)
impulse(sys_cl,'y')
title('Close-loop system with observer')

% Step response
figure;
title('Step response');
subplot(2,2,1)
impulse(sys,'r');
title('Original system')
subplot(2,2,2)
step(sysCon,'g');
title('System with controller')
subplot(2,2,3)
step(sysObs,'b')
title('system with observer')
subplot(2,2,4)
step(sys_cl,'y')
title('Close-loop system with observer')
