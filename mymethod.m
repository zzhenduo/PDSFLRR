clc;
clear;

addpath tools
addpath H:\zzd\scRNA_data
addpath(genpath('manopt'));
dataset = {'Treutlin','Goolam','mECS','Engel','Pollen','Darmanis','islet','Kolod','Tasic','Zeisel','Buettner','Macosko'};
for num=1:1
    load(['Data_' dataset{num}]);

fea=double(in_X);
[X,] = FilterGenesZero(fea);
X = normalize(X');
gnd=true_labs(:);
K = max(gnd); 
class_labels = zeros(1, K);
lambdas = [1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1];
alphas=[1e-5,1e-4,1e-3,1e-2,1e-1,1];
betas = [1:10];
beta1 = 1;%0,ablation
XX = X;
acc_index = 1;
nmi_index = 1;
ari_index = 1;
                    acc_max = 0;
                    nmi_max = 0;
                    ari_max = -1;
    for P_i = 0.1 % 0.5 
        for lambda_2 = 1.5  % 1.1:0.1:2.5,d
            for lmd_idx = 1 : length(lambdas)
                lambda = lambdas(lmd_idx); 
                for alpha_idx=1:length(alphas)
                alpha=alphas(alpha_idx);
                tic;
                Z = PDSFLRR(normc(XX),lambda,lambda_2,alpha,P_i,beta1);

                for beta_idx = 1 : length(betas)
                    beta = betas(beta_idx);
                    [U, s, ~] = svd(Z, 'econ');
                    s = diag(s);
                    r = sum(s>1e-6);
                    U = U(:, 1 : r);
                    s = diag(s(1 : r));
                    M = U * s.^(1/2);
                    mm = normr(M);
                    rs = mm * mm';
                    L = rs.^(2 * beta);
                    time_cost = toc;
                    actual_ids = spectral_clustering(L, K);
                    acc(acc_index) = accuracy(gnd, actual_ids);
                    
                    nmi(nmi_index) = Cal_NMI(gnd,actual_ids);
                   
                    ari(nmi_index) = Cal_ARI(gnd,actual_ids);
                   
                    fprintf("lambda_2= %f,lambda=%f,alpha=%f,beta=%d,NMI=%f,ARI=%f,ACC=%f\n",lambda_2,lambda,alpha,beta,nmi(nmi_index),ari(ari_index),acc(acc_index));

                    if(max(acc)>acc_max)
                        acc_max = max(acc);
                    end
                    if(max(nmi)>nmi_max)
                        nmi_max = max(nmi);
                    end
                    if(max(ari)>ari_max)
                        ari_max = max(ari);
                    end
                    acc_index = acc_index+1;
                    nmi_index = nmi_index+1;
                    ari_index = ari_index+1;
                end
                end
end
    end
    end
data = [acc(:), nmi(:), ari(:)]; 
max_values=zeros(42,3);
for col=1:3
    for group=1:42
        start_idx = (group - 1) * 10 + 1; 
        end_idx = group * 10; 
        max_values(group, col) = max(data(start_idx:end_idx, col)); 
    end
end 
file_name = 'C:/Users/ZZD/Desktop/parasPDSFLRR/Treutlin.xlsx';
writematrix(max_values, file_name);
end
% end
% plot(acc)
