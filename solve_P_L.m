% % This function was created by PDRLRR：
% % Chen H, Chen X, Tao H, et al. PDRLRR: A novel low-rank representation with projection distance regularization via manifold optimization for clustering[J]. Pattern Recognition, 2024, 149: 110198.
function [W] = solve_P_L(W, X, Z, E, Y1, mu, P_i, Lz, beta)
XE = X-X*Z;
% grad = X1*Y1' + mu*X1*(P1*X1-E)';
[d k] = size(W);
% Cxx = X*X' + P_i*eye(d);
Cxx = X*X' + P_i*eye(d);
D = W'*X;
%%
    problem.M  = stiefelgeneralizedfactory(d, k,Cxx); 
    % problem.M  = grassmanngeneralizedfactory(n, p,Cxx);
    problem.cost = @cost;
    function f = cost(W)
%     f = trace((W'*XE-E)*Y1') + mu*norm((W'*XE-E),'fro')/2;
        f = beta*trace(D*Lz*D')+trace((W'*XE-E)*Y1') + mu*norm((W'*XE-E),'fro')/2;
    end
    problem.egrad = @egrad;
    function G = egrad(W)
        G = XE*Y1' + mu*XE*(W'*XE-E)' + 2*beta*D*Lz*D';
    end
    [W, costw, info1, options1] = conjugategradient(problem, W);
%     W = problem.M.proj(W, W);
% P = W;
end

