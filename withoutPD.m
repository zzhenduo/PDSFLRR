clc;
clear;

addpath tools
% addpath history
addpath E:\ZZD\scRNAData
addpath E:\ZZD\recode\LRRNC\code
addpath(genpath('manopt'));
dataset = {'Ramskold','Treutlin','Ting','Goolam','Deng','mECS','Engel','Pollen','Darmanis','islet','Usoskin','Kolod','Tasic','Zeisel','Buettner','Macosko','Bre','Ginhoux','Grover','Leng'};
for num=14:14
    load(['Data_' dataset{num}]);

fea=double(in_X);
% X = fea';
[X,] = FilterGenesZero(fea);
X = normalize(X');
gnd=true_labs(:);
K = max(gnd); 
class_labels = zeros(1, K);
for idx =  1 : K
    class_labels(idx) = length(find(gnd == idx));
    % 统计每个类别有多少
end
lambdas = [1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1];%:0.1:1;%[0.9];         % 0.9 
betas = [1:10];
beta1 = 0;
alphas=[1e-5,1e-4,1e-3,1e-2,1e-1,1];
XX = X;
acc_index = 1;
nmi_index = 1;
ari_index = 1;
        % result_index=1;
% result= zeros(3,25);
                    acc_max = 0;
                    nmi_max = 0;
                    ari_max = 0;
% for p = 1:0.5:10%   1  10               % 10
    for P_i = 0.5%:0.01:0.1       % 0.01 
        for lambda_2 = 1.5
            % 1.94    
            for lmd_idx = 5:length(lambdas)
                lambda = lambdas(lmd_idx); 
                for alpha_idx=1:length(alphas)
                alpha=alphas(alpha_idx);
                tic;
                % [Z, E, Y1, P, dim] = LRRARD_L_Sp(normc(XX), lambda, P_i, lambda_2, beta1, p);
                % Z = LRRARD_L_Sp_ablation(normc(XX), lambda, P_i, lambda_2, beta1, p);
                Z = PDSFLRR_PD(normc(XX),lambda,lambda_2,alpha,P_i,beta1);
                % L=(Z+Z')/2;
                % % beta=0;

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
                   
                    fprintf("lambda_2=%f,lambda= %f,alpha=%f,beta=%d,NMI=%f,ARI=%f,ACC=%f\n",lambda_2,lambda,alpha,beta,nmi(nmi_index),ari(ari_index),acc(acc_index));
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
    % result;
% end
acc_max 
nmi_max 
ari_max 
end
% plot(acc)
