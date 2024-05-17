 % Written by Rupak Datta Aug 2020 %%%%%%%%%
clear, clc; 
disp('******************* LMI *******************')
warning('off','YALMIP:strict') ;
%%

n=2;

%% Comparison-1

A=[-2     0 ;                   
   0    -0.9];

B=[-1   0 ;                   
  -1    -1];

%%  Comparison-2

% A=[0    1;
%    1 0];
% 
% Ad=[0 1]'; K=-[2  2];
% 
% 
% B=Ad*K;

%% Comparison-3

% A=[-1   -1;
%    -1  -2];
% 
% B=[0   0 ;                   
%   -1   -2];



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tau1=10^(-5);
tau2=tau1+0.0001;
incr=1; pres=1; dp=0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Simulation parameters
Zn = zeros(n); ln = eye(n); Z2n = zeros(2*n);Z3n = zeros(3*n);Z4n = zeros(4*n);

[n n] = size(A);


k=6;
for l=1:k
e{l}=[zeros(n,(l-1)*n),eye(n), zeros(n,(k-l)*n)];
 end
zeroN=e{1}-e{1};

while min(pres)>0,

%% LMI variable & terms
P = sdpvar(n,n,'symmetric');
S1 = sdpvar(2*n,2*n,'symmetric');  %%Z1
S2 = sdpvar(2*n,2*n,'symmetric');  %%Z2

Z = sdpvar(3*n,3*n,'symmetric');   %S

Q1 = sdpvar(n,n,'symmetric');
Q2 = sdpvar(n,n,'symmetric');

H1 = sdpvar(n,3*n,'full');
H2 = sdpvar(n,3*n,'full');



X1 = sdpvar(6*n,n,'full');
X2 = sdpvar(6*n,n,'full');
X3 = sdpvar(6*n,n,'full');

Y1 = sdpvar(6*n,n,'full');
Y2 = sdpvar(6*n,n,'full');
Y3 = sdpvar(6*n,n,'full');



M1 = sdpvar(n,n,'full');
M2 = sdpvar(n,n,'full');
M3 = sdpvar(n,n,'full');

%%
U1 = sdpvar(n,n,'symmetric'); %R1
U2 = sdpvar(n,n,'symmetric');  %R2
U3 = sdpvar(n,n,'full');      %R3
U4 = sdpvar(n,n,'full');      %R4
U5 = sdpvar(n,n,'symmetric'); %R5


p=1;q=1;
 
R1=[(p*(U1)+p*(U1)'-q*(U2)-q*(U2)')                    2*q*(U2)                                      U3;
       (2*q*(U2))'                                 (-p*(U1)-p*(U1)'-q*(U2)-q*(U2)')                 U4;
         U3'                                            U4'                                        U5+U5'];

%% Definiton of Matrices  
%%%%%%%%%%% \dot(V)(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gamma1=[e{3};e{4};e{5}+e{6}];
Gamma2=[e{1};e{3};e{5}];
Gamma3=[e{1};e{4};e{6}];
Gamma4=[e{2};zeroN;e{1}];
Gamma5=[e{2};zeroN;-e{1}];
Gamma6=[e{2};e{1}];


Xi=[e{1};e{2};e{3};e{4};e{5};e{6}];
%% Zero Equation
Gamma=(e{1}'*M1+e{2}'*M2+e{3}'*M3);
Delta=(A*e{1}+B*e{3}-e{2});

%%

%%%%%%%%%%%%% Matrices definition %%%%%%%%%%%%%%%%%%%%%%%%
   
Theta1=e{1}'*P*e{2}+(e{1}'*P*e{2})'+(e{3}-e{1})'*H1*Gamma1+(+(e{3}-e{1})'*H1*Gamma1)'+(e{1}-e{4})'*H2*Gamma1+((e{1}-e{4})'*H2*Gamma1)'...
               +Xi'*X1*(e{1}-e{3})+(Xi'*X1*(e{1}-e{3}))'+Xi'*X2*e{5}+(Xi'*X2*e{5})'-2*Xi'*X3*e{5}-(2*Xi'*X3*e{5})'...
               +Xi'*Y1*(e{4}-e{1})+(Xi'*Y1*(e{4}-e{1}))'+Xi'*Y2*e{6}+(Xi'*Y2*e{6})'-2*Xi'*Y3*e{6}-(2*Xi'*Y3*e{6})'...
               +(e{3}-e{1})'*Q1*(e{1}-e{3})+(e{1}-e{4})'*Q2*(e{1}-e{4})-Gamma2'*R1*Gamma2...
               +(Gamma*Delta)+(Gamma*Delta)';

%%
% (t-t_k)
Theta2=e{2}'*Q2*(e{1}-e{4})+(e{2}'*Q2*(e{1}-e{4}))'...
              +e{2}'*H2*Gamma1+(e{2}'*H2*Gamma1)'...
              +Xi'*X3*(e{1}+e{3})+(Xi'*X3*(e{1}+e{3}))'...
              +Gamma6'*S2*Gamma6-Gamma1'*Z*Gamma1; %

%%

% (t_k+1-t)
% 
Theta3=e{2}'*Q1*(e{1}-e{3})+(e{2}'*Q1*(e{1}-e{3}))'...
              +e{2}'*H1*Gamma1+(e{2}'*H1*Gamma1)'+Gamma2'*R1*Gamma4+(Gamma2'*R1*Gamma4)'...
              +Xi'*Y3*(e{1}+e{4})+(Xi'*Y3*(e{1}+e{4}))'...
              +Gamma6'*S1*Gamma6+Gamma1'*Z*Gamma1; %

    %%  
      
    Z5n = zeros(6*n,n); Z6n = zeros(2*n,2*n); Z7n=zeros(2*n,1);
    
  % (t-tk)
      
LMI1=[Theta1+((tau2)*Theta2),                      sqrt(tau2)*[X1 X2],       tau2*sqrt(tau2)*[X3 Z5n];
       (sqrt(tau2)*[X1 X2])',                               -S1,                       Z6n;
       (tau2*sqrt(tau2)*[X3 Z5n])',                        Z6n',                      -3*S1  ];
    
%(tk+1)

LMI2=[Theta1+((tau2)*Theta3),                    sqrt(tau2)*[Y1 Y2],       tau2*sqrt(tau2)*[Y3 Z5n];
       (sqrt(tau2)*[Y1 Y2])',                               -S2,                       Z6n;
       (tau2*sqrt(tau2)*[Y3 Z5n])',                        Z6n',                      -3*S2];

lmi1 = LMI1; 
lmi2 = LMI2;

% 
Cst1 = lmi1;
Cst2 = lmi2;

%%

% Some basic conditions

    ineqs =[P>=0,S1>=0,S2>=0,Z>=0,Cst1<=0,Cst2<=0];

  
opstion = sdpsettings('solver','sedumi','verbose',0);
solution = solvesdp(ineqs,[], opstion);
   
[pres, dres] = checkset(ineqs);

    if incr==dp, break;
    else
        if min(pres)>0, P1s=double(P);
        else tau2=tau2-incr; incr=incr/10; pres=1;
        end
    end
    
    tau2=tau2+incr
    
end


fprintf('\n------------------------------------------------------------\n');
fprintf(' Result:   tau2= %g \n', tau2-incr);
fprintf('------------------------------------------------------------\n');





