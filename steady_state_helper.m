function l=steady_state_helper(l0,psi,c_l,gamma,etac,etal,w)
options=optimset('Display','off','TolX',1e-10);
l=fsolve(@(l) psi*(1-l)^(-etal)-gamma*c_l^(-etac)*l^(-etac) * w,l0,options);
end