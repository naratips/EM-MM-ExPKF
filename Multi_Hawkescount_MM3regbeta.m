function [mu_out,alpha_out,gamma_out,LL_out] = Multi_Hawkescount_MM3regbeta(Mu,Alpha,Gamma, Maxiter, dN, inode)
%---MM algorithm for M-dim Hakwes
%---gamma (decay rate) is unknown gamma_{ij} can be different------
%--dN is the format col=node , row=time
%--dN = [dN0, dN1, dN2,..., dNN]
Mu0= mean(dN(inode,:))/2;
%Gamma0 = mean(Gamma);
M = numel(Alpha);
N = length(dN)-1;
a_p = Mu0*N;
b_p = 1*N;
aa_p = 10*N/4-1; %Gamma0*N*4-1;
bb_p = 50*N/4-1; %1*N*4-1;
K= N;
gg = Gamma'.^[0:N-1]; 
hh = gg.*[0:N-1];
dNs = sum(dN(:,2:end-2),2);
mu_out = zeros(1,Maxiter);
alpha_out = zeros(M,Maxiter);
gamma_out = zeros(M,Maxiter);
mu_out(1) = Mu;
alpha_out(:,1) = Alpha;
gamma_out(:,1) = Gamma;
LL_out = zeros(1,Maxiter-1);
for i = 2:Maxiter
    SS= zeros(M,K); 
    TT=zeros(M,K);
    for j = 1:M
        S = conv(gg(j,:),dN(j,1:end-1));
        S = S(1:K);
        S = Alpha(j)*S;
        SS(j,:)=S;
        T = conv(hh(j,:),dN(j,1:end-1));
        T = T(1:K);
        T = Alpha(j)*T;
        TT(j,:)=T;
    end
    LL_out(i-1) = -sum(dN(inode,2:end).*log(N*Mu+sum(SS)))+sum(N*Mu+sum(SS));%'old' para
    AA = dN(inode,2:end)./(Mu+sum(SS));
    Mu = Mu/(N+b_p)*sum(AA)+(a_p-1)/(N+b_p);
    
    for j= 1:M
        
        c  = sum(AA.*SS(j,:));  
        a  = (1+Gamma(j))*dNs(j)/Alpha(j);
        b  = dN(j,N-1);
        aa = dNs(j)*Alpha(j)/(1+Gamma(j));
        cc = sum(AA.*TT(j,:));
        
        Alpha(j) = (-b+sqrt(b^2+4*a*c))/2/a; %update alpha_{inode,i}
        ppp = [-aa,0,aa+cc,aa_p+bb_p-cc,-aa_p];
        %Gamma(j) = (-bb+sqrt(bb^2+4*aa*cc))/2/aa; %update gamma_{inode,i}
        tmp = roots(ppp);
        Gamma(j)= min(tmp(tmp>0)); %choose the smallest positive sol.
       
    end
    mu_out(i)=Mu;
    alpha_out(:,i) = Alpha(:);
    gamma_out(:,i) = Gamma(:);
    gg = Gamma'.^[0:N-1]; 
    hh = gg.*[0:N-1];
end