function [Z] = PDSFLRR_PD(X,lambda,lambda_2,alpha,P_i,beta)
[m,n]=size(X);
D = X;
DTD = D'*D;
k_old = size(D,1);
Z=zeros(n,n);
J=Z;
E=zeros(m,n);
Y1=zeros(m,n);
Y2=zeros(n,n);

rho = 1.1;
mu = 1e-2;  
tol = 1e-6;
max_mu = 1e5;
max_iter = 2e2;
% thd = 0.001;
for iter=1:max_iter
    
    %% updata D  

    
    %% updata Z
    Z = inv(eye(n)+DTD)*( DTD - D'*E + J + ( D'*Y1 - Y2 )/mu);
    %% update J
     Q = L2_distance_1(D,D);
    Jj=(mu*Z+Y2-Q)/(2+mu);
    tau=alpha/(2+mu);
    J=sign(Jj).*max(abs(Jj)-tau,0);
    %% update E
    E=(mu/(2*lambda+mu))*(D+1/mu*Y1-D*Z);

    DDZ = D - D*Z;
     leq1 = DDZ-E;
    leq2 = Z-J;
    stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
%     if display && (iter==1 || mod(iter,50)==0 || stopC<tol)    % mod(a,b) ��a/b������   ÿ����50�����һ��???
%         disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
%             ',rank=' num2str(rank(Z,1e-4*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
%     end                %rank ����Ϊ Z �д��� 1e-4*norm(Z,2) ������ֵ�ĸ���
    if stopC<tol 
        break;
    else
        Y1 = Y1 + mu*leq1;  
        Y2 = Y2 + mu*leq2;
        mu = min(max_mu,mu*rho);  
    end
end    



end

