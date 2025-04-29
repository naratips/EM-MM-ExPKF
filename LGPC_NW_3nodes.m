function [yout,lambout,xout,gout] = LGPC_NW_3nodes(para,dt,Nstep)
%-------this run LGPC model in Mohler's paper 2013------------------
%---dt must be very small to allow Non-homogeneous Poison process---

%--parameter for x(t)----------------
mu = para(1:3)';
omega1 = para(4:6)';
eta = para(7:9)';
alpha1 = para(10)';

%--parameter for g(t)-----------------
omega2 = para(11:13)';
theta_1j = para(14:16)';
theta_2j = para(17:19)';
theta_3j = para(20:22)';
theta=[theta_1j';theta_2j';theta_3j'];

%---these are constant terms used in every loop----------
%---hence it is worth to store them----------------------
a1 = 1-omega1*dt;
a2 = 1-omega2*dt;
omega1_mu_dt = omega1.*mu*dt;
alpha1_sqdt = alpha1*sqrt(dt);

%---Set up the storage-----------------------------------
yout = zeros(3,Nstep+1); % data
lambout = zeros(3,Nstep+1);
xout = zeros(3,Nstep+1);
gout = zeros(3,Nstep+1);

%---Initialisation---------------------------------------------
lambnow = exp(mu);
ynow = poissrnd(lambnow*dt);

xnow = mu;
gnow = zeros(3,1);

xout(:,1) = xnow;
gout(:,1) = gnow;
yout(:,1) = ynow;
lambout(:,1) = lambnow;

for i=1:Nstep
    
%     x = a1(1)*((1-eta(1))*xnow(1)+eta(1)*sum(xnow([2,3])))+omega1_mu_dt(1)+alpha1_sqdt*randn();
%     g = a2(1)*gnow(1)+sum(theta_1j.*ynow);     
%     xnow(1) = x;
%     gnow(1) = g;
%     lambnow(1) = exp(x)+g;
%     ynow(1) = poissrnd(lambnow(1)*dt);
%     
%     x = a1(2)*((1-eta(2))*xnow(2)+eta(2)*sum(xnow([1,3])))+omega1_mu_dt(2)+alpha1_sqdt*randn();
%     g = a2(2)*gnow(2)+sum(theta_2j.*ynow);     
%     xnow(2) = x;
%     gnow(2) = g;
%     lambnow(2) = exp(x)+g;
%     ynow(2) = poissrnd(lambnow(2)*dt);
%     
%     x = a1(3)*((1-eta(3))*xnow(3)+eta(3)*sum(xnow([1,2])))+omega1_mu_dt(3)+alpha1_sqdt*randn();
%     g = a2(3)*gnow(3)+sum(theta_3j.*ynow);     
%     xnow(3) = x;
%     gnow(3) = g;
%     lambnow(3) = exp(x)+g;
%     ynow(3) = poissrnd(lambnow(3)*dt);
%  
%     xout(:,i+1) = xnow;
%     gout(:,i+1) = gnow;
%     lambout(:,i+1) = lambnow;
%     yout(:,i+1) = ynow;

        x = a1.*((1-eta).*xnow+eta.*(sum(xnow)-xnow))+omega1_mu_dt+alpha1_sqdt*randn(3,1);
        g = a2.*gnow+(theta*ynow);
        xnow = x;
        gnow = g;
        lambnow = exp(x)+g;
        ynow(1) = poissrnd(lambnow(1)*dt);
        ynow(2) = poissrnd(lambnow(2)*dt);
        ynow(3) = poissrnd(lambnow(3)*dt);
        
        xout(:,i+1) = xnow;
        gout(:,i+1) = gnow;
        lambout(:,i+1) = lambnow;
        yout(:,i+1) = ynow;
    
end



