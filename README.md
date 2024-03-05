# fiscal-policy-shocks
I need help with my Dynare code, I am new and I really need help with the coding.
Here is my mod file, I have 13 equations, 12 of them have cleared but 1 is failing this is market clearing condition its = 1 instead of 0, Please any assistant will be of good help.



@#define LOGUTILITY = 1
var
  a     ${A}$        (long_name='technology productivity')
  g      ${G}$       (long_name='government consumption')
  x     ${X}$        (long_name='public investment')
  b      ${B}$       (long_name='government transfers')
  y     ${Y}$        (long_name='output')
  c     ${C}$        (long_name='consumption')
  k     ${K}$        (long_name='capital')
  l     ${L}$        (long_name='labor')
  r     ${R}$        (long_name='interest Rate')
  w     ${W}$        (long_name='wage')
  iv    ${I}$        (long_name='investment')
  mc    ${MC}$       (long_name='marginal Costs')
  tau   ${TAU}$       (long_name='tax')
 
  
;

model_local_variable
  uc    ${U_t^C}$
  ucp   ${E_t U_{t+1}^C+g}$
  ul    ${U_t^L}$
  fk    ${f_t^F}$
  fl    ${f_t^L}$
;

varexo
 eps_a  
 eps_x    
 eps_g 
 eps_b 
;

parameters
  KAPPA   ${\kappa}$ (long_name='Private Discount Factor')
  OMEGA  ${\omega}$ (long_name='Private Discount Factor')
  BETA  ${\beta}$  (long_name='Discount Factor')
  DELTA ${\delta}$ (long_name='Depreciation Rate')
  DELTA_p ${\delta_p}$ (long_name='Depreciation Rate')
  GAMMA ${\gamma}$ (long_name='Consumption Utility ')
  CHI   ${\chi}$   (long_name='Labor Disutility ')
  @#if LOGUTILITY != 1
  ETAC  ${\eta^C}$ (long_name='consumption substitubility =gamma')
  ETAL  ${\eta^L}$ (long_name='Inverse Frisch Elasticity 1+phi')
  @#endif
  ALPHA 
  RHOA  
  RHOG  
  RHOX  
  RHOB 
;

% Parameter calibration

ALPHA = 0.05;
BETA  = 0.85;
DELTA = 0.22;
GAMMA = 0.2;
CHI   = 0.2;
KAPPA = 25.0;
DELTA_p = 0.22;
OMEGA = 0.3;
RHOA  = 0.43;
RHOX  = 0.36;
RHOG  =0.17;
RHOB  =0.3;


@#if LOGUTILITY == 0
ETAC  = 1;
ETAL  = 2;
@#endif


model;
%marginal utility of consumption and labor
@#if LOGUTILITY == 1
  #uc  = g*GAMMA*c^(-1);
  #ucp  = g*GAMMA*c(+1)^(-1);
  #ul = -CHI*(1-l)^(-1);
@#else
  #uc  = g*GAMMA*c^(-ETAC);
  #ucp  = g*GAMMA*c(+1)^(-ETAC);
  #ul = -CHI*(1-l)^(-ETAL);
@#endif

%marginal products of production
#ff = ALPHA*y/k(-1);
#fl = (1-ALPHA)*y/l;

[name='intertemporal optimality (Euler)'] %equation1
uc = BETA*ucp*(1-DELTA+r(+1)); 
[name='labor supply'] %equation2
w = -ul/uc;
[name='capital accumulation']%equation3
k = (1-DELTA)*k(-1) + (1-(KAPPA/2)*(iv/(iv(-1))-1)^2)*iv;
[name='market clearing']%equation5
y <= c + iv + x + g;
[name='production function'] %equation6
y = x*a*k(-1)^ALPHA*l^(1-ALPHA);
[name='marginal costs'] %equation7
mc = 1;
[name='labor demand'] %equation8
w = mc*fl;
[name='capital demand'] %equation9
r = mc*ff;
[name='fiscal policy'] %equation10
eps_g + eps_b + eps_x = tau;
[name='total factor productivity'] %equation11
log(a) = RHOA*log(a(-1)) + eps_a;
[name='public investment productivity']
log(x) = RHOX*log(x(-1)) + eps_x;
[name='government consumption productivity']
log(g) = RHOG*log(g(-1)) + eps_g;
[name='government transfer']
log(b) = RHOB*log(b(-1)) + eps_b;

end;



% ------------------------ %
% Steady State Computation %
% ------------------------ %
@#define AnalyticalSteadyState =1

@#if AnalyticalSteadyState == 1

steady_state_model;
x=1;
g=  1;
a = 1;
mc = 1;
b=1;
r = 1/BETA + DELTA -1;
F_L = (mc*ALPHA*a/r)^(1/(1-ALPHA));
w = mc*(1-ALPHA)*a*F_L^ALPHA;
IV_L = DELTA*F_L;
Y_L = a*(F_L)^ALPHA;  
C_L = Y_L - IV_L;
@#if LOGUTILITY==1
  l = GAMMA/CHI*C_L^(-1)*w/(1+GAMMA/CHI*C_L^(-1)*w);
@#else
  L0 = 1/3;
  l = rbc_steady_state_helper(L0, w,C_L,ETAC,ETAL,CHI,GAMMA);
@#endif
c  = C_L*l;
y  = Y_L*l;
iv = IV_L*l;
k  = F_L*l;
f=1;

end;

@#else


@#endif


steady;


shocks;

var eps_a= 0.3^2;
var eps_x= 0.3^2;
var eps_g= 0.3^2;
var eps_b= 0.3^2;
end;

stoch_simul;
stoch_simul(order=1, pruning, irf=20);






