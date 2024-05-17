 % Written by Rupak Datta Aug 2020 %%%%%%%%%
clear, clc; 
disp('******************* LMI *******************')
warning('off','YALMIP:strict') ;
%%

n=3;

J=0.005;
Rt=0.5;
d=1.2;
f=0.006;
Ind=0.0027;%L
rhooo=8;
lamopt=8.1;
pf=0.16;
Rs=1.13;
Cp=0.41;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


W1=0;
W2=5;

Tm=(((d)*(3.14)*(Rt*Rt*Rt*Rt*Rt)*Cp))/(2*(lamopt*lamopt*lamopt));

    
A{1}=[-(f+(Tm*(W1-W2)))/J                (1.5*rhooo*pf)/J               0;                   
      
      -(rhooo*pf)/Ind          -(Rs/Ind)                 -rhooo*(W1-W2) ;
        0                  rhooo*(W1-W2)                 -(Rs/Ind)];
    
A{2}=[-(f+(Tm*(W1+W2)))/J                (1.5*rhooo*pf)/J               0;                   
      
      -(rhooo*pf)/Ind          -(Rs/Ind)                 -rhooo*(W1+W2) ;
        0                  rhooo*(W1+W2)                 -(Rs/Ind)];
    
B{1}=[0 0; 1/Ind 0; 0 1/Ind ];   
B{2}=B{1};


C{1}=[1 0 0];
C{2}=C{1};



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%CASE-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Simulation parameters
Zn = zeros(n); Z2n = zeros(2*n);Z3n = zeros(3*n);Z4n = zeros(4*n);

[n n] = size(A{2});


k=6;
for l=1:k
e{l}=[zeros(2*n,(l-1)*2*n),eye(2*n), zeros(2*n,(k-l)*2*n)];
 end
zeroN=e{1}-e{1};



%% LMI variable & terms
P{1} = sdpvar(2*n,2*n,'symmetric');
P{2} = sdpvar(2*n,2*n,'symmetric');


S1 = sdpvar(4*n,4*n,'symmetric');  %%Z1
S2 = sdpvar(4*n,4*n,'symmetric');  %%Z2

Z= sdpvar(6*n,6*n,'symmetric');    %%S
Q1 = sdpvar(2*n,2*n,'symmetric');
Q2 = sdpvar(2*n,2*n,'symmetric');

H1 = sdpvar(2*n,6*n,'full');
H2 = sdpvar(2*n,6*n,'full');

X1 = sdpvar(6*n,2*n,'full');
X2 = sdpvar(6*n,2*n,'full');
X3 = sdpvar(2*n,2*n,'full');

Y1 = sdpvar(8*n,2*n,'full');
Y2 = sdpvar(8*n,2*n,'full');
Y3 = sdpvar(4*n,2*n,'full');

E = sdpvar(n,n,'symmetric');
K{1}=sdpvar(2,n,'full');
K{2}=sdpvar(2,n,'full');

%% OBSERVER

L{1}=sdpvar(n,1,'full');
L{2}=sdpvar(n,1,'full');

Psi{1} = sdpvar((24*n),(24*n),'symmetric');
Psi{2} = sdpvar((24*n),(24*n),'symmetric');

%%

U1 = sdpvar(2*n,2*n,'symmetric');
U2 = sdpvar(2*n,2*n,'symmetric');
U3 = sdpvar(2*n,2*n,'full');
U4 = sdpvar(2*n,2*n,'full');
U5 = sdpvar(2*n,2*n,'symmetric');

p=1;q=1;
 
R1=[(p*(U1)+p*(U1)'-q*(U2)-q*(U2)')                    2*q*(U2)                                      U3;
       (2*q*(U2))'                                 (-p*(U1)-p*(U1)'-q*(U2)-q*(U2)')                 U4;
         U3'                                            U4'                                        U5+U5'];


%%


U = sdpvar(2*n,2*n,'symmetric');
V{1} = sdpvar(2*n,2*n,'symmetric');
V{2} = sdpvar(2*n,2*n,'symmetric');


%%


tau2=0.454;

hbar{1}=0.02;
hbar{2}=0.5;
rho1=0.002;
rho2=0.1;
omega{1}=0.05;  %\chi_1
omega{2}=0.02;  %\chi_2


%% Definiton of Matrices  
%%%%%%%%%%% \dot(V)(t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




Gamma1=[e{3};e{4};e{5}+e{6}];
Gamma2=[e{1};e{3};e{5}];
Gamma3=[e{1};e{4};e{6}];
Gamma4=[e{2};zeroN;e{1}];
Gamma5=[e{2};zeroN;-e{1}];
Gamma6=[e{2};e{1}];

%%



Xi=[e{1};e{2};e{3};e{4};e{5};e{6}];

Xi_11=[e{2};e{3};e{4}];
Xi_12=e{1};
Xi_21=[e{1};e{2};e{4};e{6}];
Xi_22=[e{3};e{5}];





%%

Ze1=zeros(n,n);  Ze2=zeros(2,n); Ze3=zeros(1,n); Ze4=zeros(2,n);



Aa{1,1}=[A{1}               Ze1;
      A{1}-A{1}         A{1}]; 
  
Aa{2,2}=[A{2}               Ze1;
      A{2}-A{2}         A{2}]; 
  
Aa{1,2}=[A{1}             Ze1;
      A{1}-A{2}         A{2}];  
  
Aa{2,1}=[A{2}               Ze1;
      A{2}-A{1}         A{1}]; 

Ba{1,1}=[B{1}                          -B{1};
         ((B{1}-B{1}))                -((B{1}-B{1}))];  
     
     
Ba{2,2}=[B{2}                           -B{2};
         ((B{2}-B{2}))               -((B{2}-B{2})) ]; 
     
Ba{2,1}=[B{2}                           -B{2};
         ((B{2}-B{1}))                -((B{2}-B{1})) ]; 
     
Ba{1,2}=[B{1}                     -B{1};
         ((B{1}-B{2}))              -((B{1}-B{2})) ]; 
     

     
%      
La{1,1}=[L{1}*(C{1}-C{1}); -L{1}*C{1}]; %H_11
La{2,2}=[L{2}*(C{2}-C{2}); -L{2}*C{2}]; %H_22
La{1,2}=[L{2}*(C{2}-C{1}); -L{2}*C{2}]; %H_12
La{2,1}=[L{1}*(C{1}-C{2}); -L{1}*C{1}]; %H_21

Ca{1}=[C{1} Ze3];
Ca{2}=Ca{1};

E1=[E    Ze1;
   Ze1' E];

Cont{1}=[K{1}    Ze4;
       Ze4    K{1}];

Cont{2}=[K{1}    Ze4;
       Ze4     K{2}];
%% Zero Equation


for i=1:2
for k=1:2
for j=1:2
    
ln = eye(n);    
    
    
G=[ln ln];   % 3 by 6 matrix %

G1=[ln; ln];
    
Gamma7=(e{1}'+rho1*e{2}'+rho2*e{3}');
       

Phi{i,j,k}=Gamma7*(Aa{i,k}*E1*e{1}+Ba{i,k}*Cont{j}*e{3}+La{i,k}*G*e{3}-E1*e{2})...
           +(Gamma7*(Aa{i,k}*E1*e{1}+Ba{i,k}*Cont{j}*e{3}+La{i,k}*G*e{3}-E1*e{2}))';

end
end
end

%%

%%%%%%%%%%%%% Matrices definition %%%%%%%%%%%%%%%%%%%%%%%%
 

for i=1:2
for k=1:2  
for j=1:2 
Theta{i,k,j}=e{1}'*P{j}*e{2}+(e{1}'*P{j}*e{2})'+(e{3}-e{1})'*H1*Gamma1+(+(e{3}-e{1})'*H1*Gamma1)'+(e{1}-e{4})'*H2*Gamma1+((e{1}-e{4})'*H2*Gamma1)'...
               +Xi_11'*X1*(e{1}-e{3})+(Xi_11'*X1*(e{1}-e{3}))'+Xi_11'*X2*e{5}+(Xi_11'*X2*e{5})'-2*Xi_12'*X3*e{5}-(2*Xi_12'*X3*e{5})'...
               +Xi_21'*Y1*(e{4}-e{1})+(Xi_21'*Y1*(e{4}-e{1}))'+Xi_21'*Y2*e{6}+(Xi_21'*Y2*e{6})'-2*Xi_22'*Y3*e{6}-(2*Xi_22'*Y3*e{6})'...
               +(e{3}-e{1})'*Q1*(e{1}-e{3})+(e{1}-e{4})'*Q2*(e{1}-e{4})-Gamma2'*R1*Gamma2...%-Gamma3'*R2*Gamma3...
               +0.25*hbar{j}*e{1}'*V{j}*e{1}+Phi{i,k,j};
           
           
Anbu1{i,k,i}=e{1}'*P{j}*e{2}+(e{1}'*P{j}*e{2})'+(e{3}-e{1})'*H1*Gamma1+(+(e{3}-e{1})'*H1*Gamma1)'+(e{1}-e{4})'*H2*Gamma1+((e{1}-e{4})'*H2*Gamma1)'...
               +Xi_11'*X1*(e{1}-e{3})+(Xi_11'*X1*(e{1}-e{3}))'+Xi_11'*X2*e{5}+(Xi_11'*X2*e{5})'-2*Xi_12'*X3*e{5}-(2*Xi_12'*X3*e{5})'...
               +Xi_21'*Y1*(e{4}-e{1})+(Xi_21'*Y1*(e{4}-e{1}))'+Xi_21'*Y2*e{6}+(Xi_21'*Y2*e{6})'-2*Xi_22'*Y3*e{6}-(2*Xi_22'*Y3*e{6})'...
               +(e{3}-e{1})'*Q1*(e{1}-e{3})+(e{1}-e{4})'*Q2*(e{1}-e{4})-Gamma2'*R1*Gamma2...%-Gamma3'*R2*Gamma3...
               +0.25*hbar{j}*e{1}'*V{j}*e{1}+Phi{i,k,i};

Anbu2{j,k,i}=e{1}'*P{j}*e{2}+(e{1}'*P{j}*e{2})'+(e{3}-e{1})'*H1*Gamma1+(+(e{3}-e{1})'*H1*Gamma1)'+(e{1}-e{4})'*H2*Gamma1+((e{1}-e{4})'*H2*Gamma1)'...
               +Xi_11'*X1*(e{1}-e{3})+(Xi_11'*X1*(e{1}-e{3}))'+Xi_11'*X2*e{5}+(Xi_11'*X2*e{5})'-2*Xi_12'*X3*e{5}-(2*Xi_12'*X3*e{5})'...
               +Xi_21'*Y1*(e{4}-e{1})+(Xi_21'*Y1*(e{4}-e{1}))'+Xi_21'*Y2*e{6}+(Xi_21'*Y2*e{6})'-2*Xi_22'*Y3*e{6}-(2*Xi_22'*Y3*e{6})'...
               +(e{3}-e{1})'*Q1*(e{1}-e{3})+(e{1}-e{4})'*Q2*(e{1}-e{4})-Gamma2'*R1*Gamma2...%-Gamma3'*R2*Gamma3...
               +0.25*hbar{j}*e{1}'*V{j}*e{1}+Phi{j,k,i};

%%
% (t-t_k)


Theta11=e{2}'*Q2*(e{1}-e{4})+(e{2}'*Q2*(e{1}-e{4}))'...
              +e{2}'*H2*Gamma1+(e{2}'*H2*Gamma1)'...%-Gamma3'*R2*Gamma5-(Gamma3'*R2*Gamma5)'...
              +Xi_12'*X3*(e{1}+e{3})+(Xi_12'*X3*(e{1}+e{3}))'...
              +Gamma6'*S2*Gamma6-Gamma1'*Z*Gamma1; %



%%

% (t_k+1-t)
% 
Theta12=e{2}'*Q1*(e{1}-e{3})+(e{2}'*Q1*(e{1}-e{3}))'...
              +e{2}'*H1*Gamma1+(e{2}'*H1*Gamma1)'...
              +Gamma2'*R1*Gamma4+(Gamma2'*R1*Gamma4)'...
              +Xi_22'*Y3*(e{1}+e{4})+(Xi_22'*Y3*(e{1}+e{4}))'...
              +Gamma6'*S1*Gamma6+Gamma1'*Z*Gamma1; %
          

    %%  
    
    
    Z5n = zeros(2*n,2*n); Z6n = zeros(4*n,4*n); Z7n=zeros(4*n,1);  Z8n = zeros(4*n,2*n);  Z9n = zeros(2*n,2*n);
    
   
    
  % (t-tk)
      
Pratap1{i,k,j}=[Theta{i,k,j}+((tau2)*Theta11),   Xi_11'*sqrt(tau2)*[X1 X2],        Xi_12'*tau2*sqrt(tau2)*[X3 Z5n],  e{1}'*((P{1})+U),  e{1}'*((P{2})+U);
       (Xi_11'*sqrt(tau2)*[X1 X2])',                              -S1,                       Z6n,                 Z8n,                  Z8n;
       (Xi_12'*tau2*sqrt(tau2)*[X3 Z5n])',                        Z6n',                      -3*S1,                Z8n,                  Z8n;
              (e{1}'*((P{1})+U))',                                 Z8n',                      Z8n',                -V{1},                   Z9n;
              (e{1}'*((P{2})+U))',                                Z8n',                      Z8n',                Z9n',                  -V{2}];
          
Pratap2{i,k,i}=[Anbu1{i,k,i}+((tau2)*Theta11),   Xi_11'*sqrt(tau2)*[X1 X2],        Xi_12'*tau2*sqrt(tau2)*[X3 Z5n],  e{1}'*((P{1})+U),  e{1}'*((P{2})+U);
       (Xi_11'*sqrt(tau2)*[X1 X2])',                              -S1,                       Z6n,                 Z8n,                  Z8n;
       (Xi_12'*tau2*sqrt(tau2)*[X3 Z5n])',                        Z6n',                      -3*S1,                Z8n,                  Z8n;
              (e{1}'*((P{1})+U))',                                Z8n',                      Z8n',                -V{1},                   Z9n;
              (e{1}'*((P{2})+U))',                                Z8n',                      Z8n',                Z9n',                  -V{2}];
   
Pratap3{j,k,i}=[Anbu2{j,k,i}+((tau2)*Theta11),   Xi_11'*sqrt(tau2)*[X1 X2],        Xi_12'*tau2*sqrt(tau2)*[X3 Z5n],  e{1}'*((P{1})+U),  e{1}'*((P{2})+U);
       (Xi_11'*sqrt(tau2)*[X1 X2])',                              -S1,                       Z6n,                 Z8n,                  Z8n;
       (Xi_12'*tau2*sqrt(tau2)*[X3 Z5n])',                        Z6n',                      -3*S1,                Z8n,                  Z8n;
              (e{1}'*((P{1})+U))',                                 Z8n',                      Z8n',                -V{1},                   Z9n;
              (e{1}'*((P{2})+U))',                                 Z8n',                      Z8n',                Z9n',                  -V{2}];
  
%%

   
%(tk+1)


Anbalagan1{i,k,j}=[Theta{i,k,j}+((tau2)*Theta12),   Xi_21'*sqrt(tau2)*[Y1 Y2],        Xi_22'*tau2*sqrt(tau2)*[Y3 Z8n],  e{1}'*((P{1})+U),  e{1}'*((P{2})+U);
       (Xi_21'*sqrt(tau2)*[Y1 Y2])',                              -S2,                       Z6n,                 Z8n,                  Z8n;
       (Xi_22'*tau2*sqrt(tau2)*[Y3 Z8n])',                        Z6n',                      -3*S2,                Z8n,                  Z8n;
              (e{1}'*((P{1})+U))',                               Z8n',                      Z8n',                -V{1},                   Z9n;
              (e{1}'*((P{2})+U))',                              Z8n',                      Z8n',                Z9n',                  -V{2}];
          
          
Anbalagan2{i,k,i}=[Anbu1{i,k,i}+((tau2)*Theta12),   Xi_21'*sqrt(tau2)*[Y1 Y2],        Xi_22'*tau2*sqrt(tau2)*[Y3 Z8n],  e{1}'*((P{1})+U),  e{1}'*((P{2})+U);
       (Xi_21'*sqrt(tau2)*[Y1 Y2])',                              -S2,                       Z6n,                 Z8n,                  Z8n;
       (Xi_22'*tau2*sqrt(tau2)*[Y3 Z8n])',                               Z6n',                      -3*S2,                Z8n,                  Z8n;
              (e{1}'*((P{1})+U))',                         Z8n',                      Z8n',                -V{1},                   Z9n;
              (e{1}'*((P{2})+U))',                         Z8n',                      Z8n',                Z9n',                  -V{2}];


Anbalagan3{j,k,i}=[Anbu2{j,k,i}+((tau2)*Theta12),   Xi_21'*sqrt(tau2)*[Y1 Y2],        Xi_22'*tau2*sqrt(tau2)*[Y3 Z8n],  e{1}'*((P{1})+U),  e{1}'*((P{2})+U);
       (Xi_21'*sqrt(tau2)*[Y1 Y2])',                              -S2,                       Z6n,                 Z8n,                  Z8n;
       (Xi_22'*tau2*sqrt(tau2)*[Y3 Z8n])',                        Z6n',                      -3*S2,                Z8n,                  Z8n;
              (e{1}'*((P{1})+U))',                                  Z8n',                      Z8n',                -V{1},                   Z9n;
              (e{1}'*((P{2})+U))',                                  Z8n',                      Z8n',                Z9n',                  -V{2}];
   
       
   
%% Results



   
 LMI1{i,k,j}=Pratap1{i,k,j}-Psi{i};
 LMI2{i,k,i}=( omega{i}*(Pratap2{i,k,i}))-(omega{i}*Psi{i})+Psi{i};
 LMI3{i,k,j}=(omega{j}*Pratap1{i,k,j})+(omega{i}*Pratap3{j,k,i})...
                  -(omega{j}*Psi{i})-(omega{i}*Psi{j})...
                  +(Psi{i}+Psi{j});  
              


 LMI4{i,k,j}=(Anbalagan1{i,k,j})-Psi{i};
LMI5{i,k,i}=(omega{i}*(Anbalagan2{i,k,i}))-(omega{i}*Psi{i})+Psi{i};
LMI6{i,k,j}=(omega{j}*Anbalagan1{i,k,j})+(omega{i}*Anbalagan3{j,k,i})...
                  -(omega{j}*Psi{i})-(omega{i}*Psi{j})...
                  +(Psi{i}+Psi{j}); 
%  
end
end
end

%%
% 
lmi1 = LMI1{1,1,1};
lmi2 = LMI1{1,1,2}; 
lmi3 = LMI1{1,2,1};
lmi4 = LMI1{2,1,1}; 
lmi5 = LMI1{2,2,1};
lmi6 = LMI1{2,1,2};
lmi7 = LMI1{1,2,2};
lmi8 = LMI1{2,2,2};

lmi9 = LMI2{1,1,1};
lmi10 = LMI2{1,1,2}; 
lmi11 = LMI2{1,2,1};
lmi12 = LMI2{2,1,1}; 
lmi13 = LMI2{2,2,1};
lmi14 = LMI2{2,1,2};
lmi15 = LMI2{1,2,2};
lmi16 = LMI2{2,2,2};

lmi17 = LMI3{1,1,1};
lmi18 = LMI3{1,1,2}; 
lmi19 = LMI3{1,2,1};
lmi20 = LMI3{2,1,1}; 
lmi21 = LMI3{2,2,1};
lmi22 = LMI3{2,1,2};
lmi23 = LMI3{1,2,2};
lmi24 = LMI3{2,2,2};

lmi25 = LMI4{1,1,1};
lmi26 = LMI4{1,1,2}; 
lmi27 = LMI4{1,2,1};
lmi28 = LMI4{2,1,1}; 
lmi29 = LMI4{2,2,1};
lmi30 = LMI4{2,1,2};
lmi31 = LMI4{1,2,2};
lmi32 = LMI4{2,2,2};

lmi33 = LMI5{1,1,1};
lmi34 = LMI5{1,1,2}; 
lmi35 = LMI5{1,2,1};
lmi36 = LMI5{2,1,1}; 
lmi37 = LMI5{2,2,1};
lmi38 = LMI5{2,1,2};
lmi39 = LMI5{1,2,2};
lmi40 = LMI5{2,2,2};

lmi41 = LMI6{1,1,1};
lmi42 = LMI6{1,1,2}; 
lmi43 = LMI6{1,2,1};
lmi44 = LMI6{2,1,1}; 
lmi45 = LMI6{2,2,1};
lmi46 = LMI6{2,1,2};
lmi47 = LMI6{1,2,2};
lmi48 = LMI6{2,2,2};

% 
Cst1 = lmi1;
Cst2 = lmi2;
Cst3 = lmi3;
Cst4 = lmi4;
Cst5 = lmi5;
Cst6 = lmi6;
Cst7 = lmi7;
Cst8 = lmi8;

Cst9 = lmi9;
Cst10 = lmi10;
Cst11 = lmi11;
Cst12 = lmi12;
Cst13 = lmi13;
Cst14 = lmi14;
Cst15 = lmi15;
Cst16 = lmi16;

Cst17 = lmi17;
Cst18 = lmi18;
Cst19 = lmi19;
Cst20 = lmi20;
Cst21 = lmi21;
Cst22 = lmi22;
Cst23 = lmi23;
Cst24 = lmi24;

Cst25 = lmi25;
Cst26 = lmi26;
Cst27 = lmi27;
Cst28 = lmi28;
Cst29 = lmi29;
Cst30 = lmi30;
Cst31 = lmi31;
Cst32 = lmi32;

Cst33 = lmi33;
Cst34 = lmi34;
Cst35 = lmi35;
Cst36 = lmi36;
Cst37 = lmi37;
Cst38 = lmi38;
Cst39 = lmi39;
Cst40 = lmi40;

Cst41 = lmi41;
Cst42 = lmi42;
Cst43 = lmi43;
Cst44 = lmi44;
Cst45 = lmi45;
Cst46 = lmi46;
Cst47 = lmi47;
Cst48 = lmi48;


%%

% Some basic conditions
  
ineqs =[P{1}>=0,P{2}>=0,S1>=0,E1>=0,S2>=0,Z>=0,V{1}>=0,V{2}>=0,Cst1<=0,Cst2<=0,Cst3<=0,Cst4<=0,Cst5<=0,Cst6<=0,Cst7<=0,Cst8<=0,...
    Cst9<=0,Cst10<=0,Cst11<=0,Cst12<=0,Cst13<=0,Cst14<=0,Cst15<=0,Cst16<=0,Cst17<=0,Cst18<=0,...
    Cst19<=0,Cst20<=0,Cst21<=0,Cst22<=0,Cst23<=0,Cst24<=0,Cst25<=0,Cst26<=0,Cst27<=0,Cst28<=0,Cst29<=0,Cst30<=0,...
    Cst31<=0,Cst32<=0,Cst33<=0,Cst34<=0,Cst35<=0,Cst36<=0,Cst37<=0,Cst38<=0,Cst39<=0,Cst40<=0,Cst41<=0,Cst42<=0,...
    Cst43<=0,Cst44<=0,Cst45<=0,Cst46<=0,Cst47<=0,Cst48<=0];

ops = sdpsettings('verbose',0,'solver','sedumi');
yalmipdiagnostics=solvesdp(ineqs,[],ops);
% yalmipdiagnostics

[m p]=checkset(ineqs);
tmin=min(m)
if m>0  
   ['Have feasible solutions']
else
   ['no feasible solution']
end


%% Output
% Control gain

s31=double(E);
s32=double(E);

K1=double( K{1}*(inv(s31)))
K2=double(K{2}*(inv(s32)))

%% SVD Technique (Observer gain)

 s41=double(L{1});
 s42=double(L{2});

[O,S,V] = svd(C{1}); 
  
  S3=S(1,1);
   
  F=inv(V)*E*inv(V)';
  F1=F(1,1);
  
  

LL1=s41*O*S3*double(F1)*inv(double(S3))*inv(O)
LL2=s42*O*S3*double(F1)*inv(double(S3))*inv(O)
  
