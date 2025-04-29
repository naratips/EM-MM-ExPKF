rng(1)
%--this version estimates only mu,omega1---------------------------
%----------(1) generate data with "true" parameters--------
%--parameter for x(t)----------------
dt=0.1;
mu = [.5 .5 .5];
omega1 = [0.25, 0.25, 0.25]/dt;
eta = [0.1, 0.1, 0.1]*3;
alpha1 = 0.125;

%--parameter for g(t)-----------------
omega2 = [0.9, 0.9, 0.9]/dt;
theta_1j = [0.75 0 1]*2;
theta_2j = [0.5 0.5 0]*2;
theta_3j = [0 0 1]*2;
dt=.1; 
truepara= [mu,omega1,eta,alpha1,omega2,theta_1j,theta_2j,theta_3j];

Nstep = 1000;
[dN,lambtrue,xtrue,gtrue] = LGPC_NW_3nodes(truepara,dt,Nstep);

% 
% parain = [mu*1.25, omega1/2,eta*0.5, alpha1, omega2/2,...
%           theta_1j*1.5, theta_2j*1.5, theta_3j*1.5]; %
parain0 = [mu*1.25,omega1*1.5,eta*0.9,alpha1,omega2*1.5,0.9*ones(1,9)];
%parain0 = truepara;
parain = parain0;

Nem=5;
options = optimoptions('fminunc','Display','none','Algorithm','quasi-newton');
options_con = optimoptions('fmincon','Display','none','Algorithm','active-set');

%options = optimoptions('fminunc','Display','none','Algorithm','trust-region','SpecifyObjectiveGradient',true);
err = norm(parain-truepara)./norm(truepara);
paras =zeros (numel(parain),Nem+1);
paras(:,1) = parain(:);
err_lamb = [];

for j=1:Nem
    j
%======E-STEP============================================================
%---Assuming that all parameters are known and initialize samples--------
Nsamp = 300;
state_ini= parain(1:3)+10*alpha1*sqrt(dt)*randn(Nsamp,3);
[xsmooth,gsmooth,lambsmooth,Wsmooth,xout,gout,lambout,Wout]=PFsmoothing_LGPC_3nodes(dN,state_ini,parain,dt);
err_lamb = [err_lamb norm(lambtrue-mean(mean(lambsmooth,3)'))/norm(lambtrue)];

%======M-STEP=============================================================
%--using smoothing particles to minimise the parameter-------------------
fx1 = @(p)obj_EM3nodes_xt_novar(p,dt,xsmooth,1);
fx2 = @(p)obj_EM3nodes_xt_novar(p,dt,xsmooth,2);
fx3 = @(p)obj_EM3nodes_xt_novar(p,dt,xsmooth,3);
fg1= @(p)obj_EM3nodes_dN(p,dt,dN,squeeze(xsmooth(:,:,1)),1);
fg2= @(p)obj_EM3nodes_dN(p,dt,dN,squeeze(xsmooth(:,:,2)),2);
fg3= @(p)obj_EM3nodes_dN(p,dt,dN,squeeze(xsmooth(:,:,3)),3);

p1 = parain([1,4,7]);
p2 = parain([2,5,8]);
p3 = parain([3,6,9]);
p4 = parain([11,14,15,16]);
p5 = parain([12,17,18,19]);
p6= parain([13,20,21,22]);
BB=zeros(3,3);
CC=zeros(4,3);
    parfor K=1:6
        if K==1
        lb = zeros(1,3); ub = [2,1/dt, 1];
        [paranewx1,~] = fmincon(fx1,p1,[],[],[],[],lb,ub,[],options_con);
        BB(:,K)=paranewx1;
        elseif K==2
        lb = zeros(1,3); ub = [2,1/dt, 1];
        [paranewx2,~] = fmincon(fx2,p2,[],[],[],[],lb,ub,[],options_con);
        BB(:,K)=paranewx2;
        elseif K==3
        lb = zeros(1,3); ub = [2,1/dt, 1];
        [paranewx3,~] = fmincon(fx3,p3,[],[],[],[],lb,ub,[],options_con);
        BB(:,K)=paranewx3;
        elseif K==4
        lb = zeros(1,4); ub = [1/dt, 5 , 5, 5];
        [paranewg1,~] = fmincon(fg1,p4,[],[],[],[],lb,ub,[],options_con);
        CC(:,K)=paranewg1(:);
        elseif K==5
        lb = zeros(1,4); ub = [1/dt, 5 , 5, 5];
        [paranewg2,~] = fmincon(fg2,p5,[],[],[],[],lb,ub,[],options_con);
        CC(:,K)=paranewg2(:);
        elseif K==6
        lb = zeros(1,4); ub = [1/dt, 5 , 5, 5];
        [paranewg3,~] = fmincon(fg3,p6,[],[],[],[],lb,ub,[],options_con);
        CC(:,K)=paranewg3(:);
        end
    end %end parfor
CC(:,1:3)=[];
parain([1,4,7]) = BB(:,1);
parain([2,5,8]) = BB(:,2);
parain([3,6,9]) = BB(:,3);
parain([11,14,15,16]) = CC(:,1);
parain([12,17,18,19]) = CC(:,2);
parain([13,20,21,22]) = CC(:,3);

err=[err norm(parain-truepara)/norm(truepara)];
paras(:,j+1)=parain(:);

end

figure
tplot0= 100; tplot1=500;
subplot(311)
kk=1;
plot(squeeze(lambsmooth(1:end,:,kk)),'r-'), hold on
plot(lambtrue(kk,:),'k--','linewidth',2), hold on
bar(dN(kk,:))
xlim([tplot0 tplot1])
xlabel('step','fontsize',14)
ylabel(['Intensity of Node ',num2str(kk)])
subplot(312)
kk=2;
plot(squeeze(lambsmooth(1:end,:,kk)),'r-'), hold on
plot(lambtrue(kk,:),'k--','linewidth',2), hold on
bar(dN(kk,:))
xlim([tplot0 tplot1])
xlabel('step','fontsize',14)
ylabel(['Intensity of Node ',num2str(kk)])
subplot(313)
kk=3;
plot(squeeze(lambsmooth(1:end,:,kk)),'r-'), hold on
plot(lambtrue(kk,:),'k--','linewidth',2), hold on
bar(dN(kk,:))
xlim([tplot0 tplot1])
xlabel('step','fontsize',14)
ylabel(['Intensity of Node ',num2str(kk)])
set(gcf,'color','w')


figure, subplot(121)
heatmap(reshape(parain(end-8:end),3,3)')
title('Approximated influence matrix')
subplot(122)
heatmap(reshape(truepara(end-8:end),3,3)')
title('Truth')
set(gcf,'color','w')



figure
subplot(311)
kk=1;
plot(squeeze(xsmooth(1:end,:,kk)),'r-'), hold on
plot(xtrue(kk,:),'k--','linewidth',2), hold on

xlim([tplot0 tplot1])
xlabel('step','fontsize',14)
ylabel(['Intensity of Node ',num2str(kk)])
subplot(312)
kk=2;
plot(squeeze(xsmooth(1:end,:,kk)),'r-'), hold on
plot(xtrue(kk,:),'k--','linewidth',2), hold on

xlim([tplot0 tplot1])
xlabel('step','fontsize',14)
ylabel(['Intensity of Node ',num2str(kk)])
subplot(313)
kk=3;
plot(squeeze(xsmooth(1:end,:,kk)),'r-'), hold on
plot(xtrue(kk,:),'k--','linewidth',2), hold on

xlim([tplot0 tplot1])
xlabel('step','fontsize',14)
ylabel(['Intensity of Node ',num2str(kk)])
set(gcf,'color','w')
