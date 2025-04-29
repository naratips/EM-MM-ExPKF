 function [xsmooth,gsmooth,lambsmooth,Wsmooth,xout,gout,lambout,Wout] =...
           PFsmoothing_LGPC_3nodes(y,state,para,dt)
%-----Assuming parameter value is fixed----------------------------------

%--parameter for x(t)----------------
mu = para(1:3);
omega1 = para(4:6);
eta = para(7:9);
alpha1 = para(10);

%--parameter for g(t)-----------------
omega2 = para(11:13);
theta_1j = para(14:16);
theta_2j = para(17:19);
theta_3j = para(20:22);

theta=[theta_1j;theta_2j;theta_3j];
%---these are constant terms used in every loop----------
%---hence it is worth to store them----------------------
a1 = 1-omega1*dt;
a2 = 1-omega2*dt;
omega1_mu_dt = omega1.*mu*dt;
alpha1_sqdt = alpha1*sqrt(dt);

%---------------------------------------------------------------------
x = state;  % each row is a sample, 
g= zeros(1,3);
lamb = exp(x)+g;


%-------set up storage to record results-------------------------------
Nstep = size(y,2);
Nsamp = size(state,1);
xout = zeros(Nstep,Nsamp,3); 
gout = zeros(Nstep,3); 
lambout = zeros(Nstep,Nsamp,3);
%----------------------------------------------------------------------
W = ones(Nsamp,1)/Nsamp; %---particle weight--------
Wout = zeros(Nstep,Nsamp);
resamp_thresh = 0.5;

%--resampling parameter---------
ss= sqrt(0.2);
Np = 1/sqrt(Nsamp);
%--(1)------Filtering step-------------------------------------------
for i=1:Nstep 
        dN_now = y(:,i);
%--log likelihood----------------------------------------------------
        LH = dN_now(:)'.*log(lamb)-lamb*dt; %--log likelihood of Poisson
        log_Wnew = log(W)+sum(LH,2);
        log_Wnew  = log_Wnew-max(log_Wnew);
        Wtmp =  exp(log_Wnew); 
        W =  Wtmp/sum(Wtmp);
        if   1/sum(W.^2)/Nsamp < resamp_thresh %test for resampling            
             sampIndex = ResampSimp(W,Nsamp);
             x = x(sampIndex,:)+ss*std(x)*Np.*randn(Nsamp,3);
             W= ones(Nsamp,1)/Nsamp;
        end
        Wout(i,:) = W(:)';
        xout(i,:,:) = x;
        gout(i,:) = g;
        lambout(i,:,:) = exp(x)+g;
        
%-------------"Forecast step"---------------------------------------
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
%     ynow(2) = poissrnd(lambnow(1)*dt);
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


        xnew = a1.*((1-eta).*x+eta.*(sum(x,2)-x))+omega1_mu_dt+alpha1_sqdt*randn(Nsamp,3);
        %xnew = a1.*x+omega1_mu_dt+alpha1_sqdt*randn(Nsamp,3);
        gnew = a2.*g+(theta*dN_now)';
        lambnew = exp(xnew)+gnew;

        x = xnew;
        g = gnew;
        lamb = lambnew;
        
end 

%%--(2)-------Smoothing ---------------------------------------------
gsmooth = g;
Nsmooth = ceil(Nsamp/4);
xsmooth = zeros(Nstep,Nsmooth,3);
lambsmooth = zeros(Nstep,Nsmooth,3);
sampIndex = ResampSimp(Wout(end,:)/sum(Wout(end,:)),Nsmooth);
lambsmooth(end,:,:) = lambout(end,sampIndex,:);
xsmooth(end,:,:) = xout(end,sampIndex,:);

for i = (Nstep-1):-1:1
    for j = 1:Nsmooth
        %---log Likeliood for transition probability
        x=squeeze(xout(i,:,:));
        xf = a1.*((1-eta).*x+eta.*(sum(x,2)-x))+omega1_mu_dt;
        LH = -log(alpha1)-0.5*((squeeze(xsmooth(i+1,j,:))'-xf).^2)./alpha1_sqdt/alpha1_sqdt;
        log_Wsmooth = log(Wout(i,:))+sum(LH,2);
        log_Wsmooth  = log_Wsmooth-max(log_Wsmooth);
        Wtmp =  exp(log_Wsmooth); 
        Wsmooth =  Wtmp/sum(Wtmp);
        sampIndex = ResampSimp(Wsmooth,1);
        xsmooth(i,j,:) = xout(i,sampIndex,:);
        lambsmooth(i,j,:) = lambout(i,sampIndex,:);
        
    end
end

%disp('done')

% 
% %---1B::--Do BSPS-------------------------------------------------------
% nsmooth=10;
% sampIndex = ResampSimp(Wp(:,end),nsmooth);
% xs = zeros(nsmooth,Nstep);
% xs(:,end)= xp(sampIndex,end);
% for j=(Nstep-1):-1:1
%     for k=1:nsmooth
%     %--log likelihood----------------------------------------------------
%     LH = -0.5*(xs(k,j+1)-phi*xp(:,j)).^2/Q;
%     log_Wnew = log(Wp(:,j))+LH;
%     log_Wnew  = log_Wnew-max(log_Wnew);
%     Wtmp =  exp(log_Wnew); 
%     W =  Wtmp/sum(Wtmp);
%     sampIndex = ResampSimp(W,1);
%     xs(k,j) = xp(sampIndex,j);
%     end
% end

