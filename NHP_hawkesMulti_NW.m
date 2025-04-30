function dNout = NHP_hawkesMulti_NW(dt,mu,alpha,omega,Nstep)
% this use Non-homogeneous version of N-dimensional Hawkes
% all parameters are assumed to be constant
% Use it to generate the truth run
% omega: in the paper it is gamma
mu=mu(:);
M = length(mu);  % no. of nodes
%-------Check stability------------------------------
% omega = 1-beta*dt;
A = alpha/((1-omega)/dt);
dd=eigs(A,1)
if dd >=1
    disp('Unstable process. Stop')
    return
end
%--------initializing the first step-----------------------
out = zeros(M,Nstep);
lamb = mu;
out(:,1) = lamb;
dNout=sparse(M,Nstep);
dN_now = poissrnd(lamb*dt);
dNout(:,1) = dN_now(:);

for i=2:Nstep
%---predicted "observation step"---------------------------------------
lamb=mu+omega*(lamb-mu)+sum(alpha.*dN_now(:)',2);
out(:,i)=lamb;
dN_now = poissrnd(lamb*dt);
dNout(dN_now>0,i) = dN_now(dN_now>0);
end


