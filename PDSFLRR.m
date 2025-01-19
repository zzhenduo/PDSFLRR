function [Z] = PDSFLRR(X,lambda,lambda_2,alpha,P_i,beta)
[m,n]=size(X);
P0 = eye(m);
D = P0'*X;
k_old = size(D,1);
Z=zeros(n,n);
J=Z;
E=zeros(m,n);
Y2=zeros(n,n);

rho = 1.1;
mu = 1e-2;  
tol = 1e-6;
max_mu = 1e5;
max_iter = 2e2;
for iter=1:max_iter
     %% updata P  
    [~, s, V] = svd(D, 'econ');
    s = diag(s);
    lmd = 1/lambda_2;
    k_new = length(find(s > lmd));%%d
    if iter <= 1 || k_old ~= k_new
        X_1 = bsxfun(@minus,X, mean(X,2));
        C = cov(X_1*X_1');
        [eigvector,~]=eigs(C,k_new); 
        X_1 = eigvector'*X_1;
        P0 = P0(1 : k_new, 1 : k_new);
        E = sparse(size(X_1,1),n);
        Y1 = zeros(size(X_1,1),n);
    end
    Wz = (J+J')/2;
    Dz = diag(sum(Wz));
    Lz = Dz-Wz;
    P = solve_P_L(P0, X_1, Z, E, Y1, mu, P_i, Lz, beta);
    P0 = P(:, 1 : k_new);
    dim(iter) = k_new;
    if k_old > k_new
        Y1 = Y1(1:k_new,:);
        E = E(1:k_new,:);
    elseif k_old < k_new
        Y1(k_old+1:k_new,:) = 0;
        E(k_old+1:k_new,:) = 0;
    end
    k_old = k_new;
    %% updata D  
    D = P0'*X_1;
    DTD = D'*D;
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
    if stopC<tol 
        break;
    else
        Y1 = Y1 + mu*leq1;  
        Y2 = Y2 + mu*leq2;
        mu = min(max_mu,mu*rho);  
    end
end    



end

