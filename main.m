clear; clc;

%% Load Your Data here %%

load('')

%% For Index-I system %%

% A= J1-J2*(J4\J3);
% E=E1;
% B= B1-J2*(J4\B2);
% C=C1-C2*(J4\J3);
%  D = - (C2*(J4\B2));



E = eye(n);
D = 0;
iter = 30; 
dim = 15; 
tol = 1e-15; tolY=1e-12;

lp =0;
hp = 5; 
tp = 200; 
s1 = eigs(-A,1,'lm');s2 = eigs(-A,1,'sm');

%% GA Parameters %%
gaParam.varmax = s1; gaParam.varmin = s2; nVar = 10;
gaParam.Varsize = [1,nVar];
gaParam.dp = 30; 
gaParam.npop = 10;
gaParam.MaxIt = 20;
gaParam.pC = 1;
gaParam.beta = 1; 
gaParam.mu = 0.001; 
gaParam.sigma = 0.001;
gaParam.gamma= 0.1;


%%%%%%%RKSM-GA%%%%%%%%%%
[gC,resnorm1,Trksm_ga_C] = rksm_ga(E,A,B,iter,n,tol,tolY,gaParam);
[gO,resnorm2,Trksm_ga_O] = rksm_ga(E',A',C',iter,n,tol,tolY,gaParam);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%ADI-GA%%%%%%%%%%
[gCa,resnorm3,Tadi_ga_C] = adi_ga(E,A,B,iter,tol,gaParam);
[gOa,resnorm4,Tadi_ga_O] = adi_ga(E',A',C',iter,tol,gaParam);
%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%BT%%%%%%%%%%
[gBal_A,gBal_B,gBal_C,gBal_E] = BT(E,A,B,C,gC,gO,dim);
[gBal_Aa,gBal_Ba,gBal_Ca,gBal_Ea] = BT(E,A,B,C,gCa,gOa,dim);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Error Calculation %%%%%%%%%%
[err] = Error(A,B,C,E,D,gBal_A,gBal_B,gBal_C,gBal_E,gBal_Aa,gBal_Ba,gBal_Ca,gBal_Ea); 
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Plot%%%%%%%%%%
[space,O_hank,ga_r_hank,ga_abs_err,ga_rel_err] =tf_plot(A,B,C,E,D,gBal_A,gBal_B,gBal_C,gBal_E,gBal_Aa,gBal_Ba,gBal_Ca,gBal_Ea,
    lp,hp,tp,resnorm1,resnorm2,resnorm3,resnorm4,resnorm5,resnorm6,resnorm7,resnorm8,resnorm9,resnorm10,resnorm11,resnorm12);
%%%%%%%%%%%%%%%%%%%%%%
