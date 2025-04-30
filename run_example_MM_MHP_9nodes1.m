% Example: Running MM algorithm for a non-homogeneous Poisson process
% modeled by a discrete-time Hawkes model with expoential decay kernel.

rng(1)

mu = [5:-0.4:4 .5:-0.04:.3];
K = length(mu);
beta = .175;

alpha =spdiags(ones(K,1)*[0.1 0.2 .4 .5 .4 0.2 0.1],-3:3,K,K)/2;
msk = zeros(K,K); msk(1:3,1:3)=1; msk(4:K,4:K)=1;
alpha = alpha.*msk;




Nstep=1000;
dt=1;
dNout = NHP_hawkesMulti_NW(dt,mu,alpha,beta,Nstep);

Ntotal = sum(dNout);


%---start the process after the first observation-----------------------
id1 = 1;
while 1
    if nnz(dNout(:,id1))==0
        id1=id1+1; 
    else
        break
    end
end
if id1>1, dNout(:,1:id1-1)=[]; end
Nstep=size(dNout,2);
%----initialisation------------------------------------------------------

Aout = zeros(9,9);

%----Run over each node---------------------------------------------------
for inode = 1:K

[mm,aa,bb,ll1]=Multi_Hawkescount_MM3regbeta(mu(inode)/4,ones(1,K),0.1*ones(1,K),500, full(dNout),inode); 

Aout(inode,:)= aa(:,end)';
end



figure(1)
cc=colormap(gray);
cc=flipud(cc);
subplot(1,2,1)
imagesc(Aout)
xticks(1:K)
yticks(1:K)
axis square
title(['Data length = ', num2str(Nstep)])
txtp= num2cell(Aout);
set(gcf,'color','w')
colormap(cc);
colorbar('southoutside')
ax = gca;
ax.FontSize = 6; 

subplot(1,2,2)
imagesc(alpha)
xticks(1:K)
yticks(1:K)
axis square
title('Ground truth')
txtp= num2cell(Aout);
set(gcf,'color','w')
colormap(cc);
colorbar('southoutside')
ax = gca;
ax.FontSize = 6; 


