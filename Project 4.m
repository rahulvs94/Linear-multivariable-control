%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Multivariable Control
%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) Design a state-feedback controller to place the eigenvalues along the real axis between
% -3 and -8 with imaginary components less than + or - j
% (b) Construct an observer (open or closed-loop) to implement the state feedback
% (c) Find the model for the entire system – including the observer and feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 
clear all
clc

A = [-3 1 -2 -1 1 0 0 -1;
    0 -5 0 0 0 0 0 0;
    1 0 -3 -2 0 0 1 -1;
    -1 2 -1 -4 -1 -1 0 0;
    1 0 1 0 -2 0 -1 2;
    0 0 -1 1 2 -5 0 -1;
    0 1 -1 0 0 0 -5 0;
    1 0 1 -1 -1 0 -1 -3];
B = [1 0;
    0 0;
    0 0;
    0 0;
    0 -1;
    -1 -1;
    0 -1;
    1 -1];
C = [1 0 1 0 1 1 0 1;
    0 1 1 1 1 0 1 0];

% State-space representation
sys = ss(A,B,C,0);

% Check controllability
P = ctrb(A,B);
uncontrollable_states = length(A)-rank(P)

% Removing uncontrollable state
sys1 = minreal(sys)
[Am,Bm,Cm,Dm] = ssdata(sys1);
 
% Pole placement 
Pc = [-3.5+i*0.5, -3.5-i*0.5, -4.5+i*0.5, -4.5-i*0.5, -5.5+i*0.5, -5.5-i*0.5, -7];
K = place(Am,Bm,Pc)
Po = 5*Pc;
L = (place(Am',Cm',Po))'

% System with state feedback controller gain
Acon = Am - Bm*K;
Bcon = Bm;
Ccon = Cm;
% State-space representation
sysCon = ss(Acon,Bcon,Ccon,0)
 
% System with observer gain
Aobs = [Am    zeros(size(Am));
        L*Cm  Am-L*Cm];
Bobs = [Bm;
        Bm];
Cobs = [Cm              zeros(size(Cm));
        zeros(size(Cm)) Cm];
% State-space representation
sysObs = ss(Aobs,Bobs,Cobs,0)
 
% Entire system with observer and state feedback controller
Acl = [Am   -Bm*K;
       L*Cm  Am-L*Cm-Bm*K];
Bcl = [Bm;
       Bm];
Ccl = [Cm zeros(size(Cm))];
% State-space representation
sysCl = ss(Acl,Bcl,Ccl,0)
 
% Step response of designed systems
figure;
subplot(1,2,1);
step(sys1,'r');
title('Minimal realized system');
subplot(1,2,2);
step(sysCl,'b');
title('Entire system');