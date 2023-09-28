@#define LOGUTILITY=0
var
y ${Y}$ (long_name='output')
c ${C}$ (long_name='consumption')
k ${K}$ (long_name='capital')
l ${l}$ (long_name='labor')
A ${A}$ (long_name='productivity')
r ${R}$ (long_name='interest rate')
w ${W}$ (long_name='wage')
iv ${I}$ (long_name='investment')
mc ${MC}$ (long_name='marginal cost')
;

model_local_variable
uc ${U_t^C}$
ul ${U_t^L}$
fk ${f_t^K}$
fl ${f_t^L}$
ucp ${Et U_{t+1}^C}$
;

varexo
epsa ${\varepsilon^A}$ (long_name='Productivity Shock')
;

parameters
beta ${\beta}$ (long_name='Discount Factor')
delta ${\delta}$ (long_name='Depreciation Rate')
gamma ${\gamma}$ (long_name='Consumption Utility Weight')
pssi ${\psi}$ (long_name='Leasure Utility Weight')
@#if LOGUTILITY==0
etac ${\eta^C}$ (long_name='Risk Aversion')
etal ${\eta^L}$ (long_name='Inverse Frisch Elasticity')
@#endif
alpha ${\alpha}$ (long_name='Output elasticity of capital')
rhoa ${\rho^A}$ (long_name='Consistency of TFP Shock')
;
%--------------------------%
% parameter calibration    %
%--------------------------%
alpha=0.35;
beta=0.99;
pssi=1.6;
delta=0.025;
gamma=1;
rhoa=0.9;
@#if LOGUTILITY==0
etac=2;
etal=1;
@#endif

%--------------------------%
%         model            %
%--------------------------%
model;
%margina utiity of consumption and labor
@#if LOGUTILITY==1
 #uc=gamma*c^(-1);
 #ucp=gamma*c(+1)^(-1);
 #ul=-pssi*(1-l)^(-1);
 @#else
 #uc=gamma*c^(-etac);
 #ucp=gamma*c(+1)^(-etac);
 #ul=-pssi*(1-l)^(-etal);
@#endif
%marginal products
#fl=(1-alpha)*y/l;
#fk=alpha*A*y/(k(-1));
[name='Euler equation']
uc=beta*ucp*(1-delta+r);
[name='labor supply']
w=-ul/uc;
[name='capital accumuation']
k=(1-delta)*k(-1)+iv;
[name='market clearing']
y=c+iv;
[name='production function']
y=A*k(-1)^alpha* l^(1-alpha);
[name='marginal costs']
mc=1;
[name='labor demand']
w=mc*fl;
[name='capital demand']
r=mc*fk;
[name='TFP stochastic process']
log(A)=rhoa*log(A(-1))+epsa;
end;

write_latex_definitions;
write_latex_parameter_table;
write_latex_original_model;
%write_latex_dynamic_model;
write_latex_static_model;
collect_latex_files;

%--------------------------%
% steady state computation %
%--------------------------%
steady_state_model;
A=1;
r=1/beta+delta-1;
mc=1;
k_l=(r/alpha)^(1/(alpha-1));
w=(1-alpha)*k_l^alpha;
iv_l=delta*k_l;
y_l=k_l^alpha;
c_l=y_l-iv_l;
@#if LOGUTILITY==1
    l=1/(1+pssi*c_l/(gamma*w));
@#else
    l0=1/3;
    l=steady_state_helper(l0,pssi,c_l,gamma,etac,etal,w);
@#endif
k=k_l*l; 
iv=iv_l*l;
y=y_l*l;
c=c_l*l;
end;

steady;



