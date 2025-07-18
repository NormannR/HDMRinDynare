addpath /home/normann/Dropbox/dynare/matlab;

dynare irbc_small_so;
%%
load("quadrature.mat");
%%
% Simulation
dr = oo_.dr;
T = 11000;
r = 100;
sim = NaN(M_.endo_nbr, T, r);
exo = NaN(M_.exo_nbr, T, r);
chol_S = chol(M_.Sigma_e);
ys = dr.ys;
y_fo = zeros(size(ys,1),T+M_.maximum_lag,r);
y_so_pruned = zeros(size(ys,1),T+M_.maximum_lag,r);
k2 = M_.nstatic+(1:M_.nspred);
order_var = dr.order_var;
order_var_k2 = order_var(k2);
ys_order_var = ys(order_var);
ys_order_var_k2 = ys(order_var_k2);
constant = ys_order_var+.5*dr.ghs2;
ghxx = dr.ghxx;
ghuu = dr.ghuu;
ghxu = dr.ghxu;
ghx = dr.ghx;
ghu = dr.ghu;
for n=1:r
    exo(:,:,n) = chol_S*randn(M_.exo_nbr,T);
    y_fo(:,1,n) = ys;
    y_so_pruned(:,1,n) = ys;
    % Pruned
    for t = 2:T+M_.maximum_lag
        yhat1 = y_fo(order_var_k2,t-1,n)-ys_order_var_k2;
        yhat2 = y_so_pruned(order_var_k2,t-1,n)-ys_order_var_k2;
        epsilon = exo(:,t-1,n);
        abcOut1 = A_times_B_kronecker_C(.5*ghxx,yhat1);
        abcOut2 = A_times_B_kronecker_C(.5*ghuu,epsilon);
        abcOut3 = A_times_B_kronecker_C(ghxu,yhat1,epsilon);
        y_so_pruned(order_var,t,n) = constant + ghx*yhat2 + ghu*epsilon ...
            + abcOut1 + abcOut2 + abcOut3;
        y_fo(order_var,t,n) = ys_order_var + ghx*yhat1 + ghu*epsilon;
    end
end
%%
tic;
% Equation errors
yhat1 = y_fo(order_var_k2,:,:) - ys_order_var_k2;
yhat2 = y_so_pruned(order_var_k2,:,:)-ys_order_var_k2;
yhat1_sq = repelem(yhat1, size(yhat1,1), 1, 1).*repmat(yhat1, size(yhat1,1), 1, 1);
eps_sq = repelem(int_nodes, size(int_nodes,1), 1, 1).*repmat(int_nodes, size(int_nodes,1), 1, 1);
yhat1_eps = permute(repmat(repelem(yhat1, size(int_nodes,1), 1, 1), 1, 1, 1, size(int_nodes,2)), [1 4 2 3]).*repmat(int_nodes, size(yhat1,1), 1, size(yhat1,2), size(yhat1,3));
abcOut1 = pagemtimes(.5*ghxx, yhat1_sq);
abcOut2 = .5*ghuu*eps_sq;
abcOut2 = repmat(abcOut2, 1, 1, size(yhat1,2), size(yhat1,3));
abcOut2 = permute(abcOut2, [1 3 4 2]);
abcOut3 = pagemtimes(ghxu, yhat1_eps);
abcOut3 = permute(abcOut3, [1 3 4 2]);
ghx_yhat2 = pagemtimes(ghx,yhat2);
ghu_eps = repmat(ghu*int_nodes, 1, 1, size(yhat1,2), size(yhat1,3));
ghu_eps = permute(ghu_eps, [1 3 4 2]);
y_next = constant + ghx_yhat2 + ghu_eps + abcOut1 + abcOut2 + abcOut3;
%%
y = NaN(3*M_.endo_nbr,1);
num_nodes = numel(int_weights);
resid = NaN(M_.endo_nbr,num_nodes);
vec_errors = NaN(M_.endo_nbr, T, r);
for t = 2:T+M_.maximum_lag
    for n=1:r
        y(1:M_.endo_nbr) = y_so_pruned(:,t-1,n);
        y(M_.endo_nbr+(1:M_.endo_nbr)) = y_so_pruned(:,t,n);
        x = exo(:,t-1,n);
        for k=1:num_nodes
            y(2*M_.endo_nbr+order_var) = y_next(:,t,n,k);
            resid(:,k) = irbc_small_so.sparse.dynamic_resid(y, x, M_.params, ys);
        end
        vec_errors(:,t-1,n) = resid*int_weights;
    end
end
toc;
%%
% In what domain do the simulations live?
% quantile(reshape(y_so_pruned,M_.endo_nbr,[]), 0.1, 2)
%%
% quantile(reshape(y_so_pruned,M_.endo_nbr,[]), 0.9, 2)
%%
% max(reshape(y_so_pruned,M_.endo_nbr,[]), [], 2)
%%
% min(reshape(y_so_pruned,M_.endo_nbr,[]), [], 2)
%%
% Get the average and quantile of EE errors
burnin=0.0909;
q=0.999;
start=round(T*burnin);
system_eqs = 1:2:M_.endo_nbr-1;
%%
errMat = abs(vec_errors(system_eqs,start:end,:));
%%
log10(mean(errMat(:)))
%%
q_err = quantile(errMat, q, [2]);
log10(mean(q_err(:)))
%%
% % Equation errors in loop form, more readable but slower than the code above
% tic;
% errors = NaN(M_.endo_nbr, T, r);
% y = NaN(3*M_.endo_nbr,1);
% x = NaN(M_.exo_nbr,1);
% y_next = NaN(M_.endo_nbr,1);
% num_nodes = numel(int_weights);
% resid = NaN(M_.endo_nbr,num_nodes);
% for t = 2:T+M_.maximum_lag
%     for n=1:r
%         y(1:M_.endo_nbr) = y_so_pruned(:,t-1,n);
%         y(M_.endo_nbr+(1:M_.endo_nbr)) = y_so_pruned(:,t,n);
%         x = exo(:,t-1,n);
%         yhat1 = y_fo(order_var_k2,t,n)-ys_order_var_k2;
%         yhat2 = y_so_pruned(order_var_k2,t,n)-ys_order_var_k2;
%         abcOut1 = A_times_B_kronecker_C(.5*ghxx,yhat1);
%         common = constant + ghx*yhat2 + abcOut1; 
%         for k=1:num_nodes
%             epsilon = int_nodes(:,k);
%             abcOut2 = A_times_B_kronecker_C(.5*ghuu,epsilon);
%             abcOut3 = A_times_B_kronecker_C(ghxu,yhat1,epsilon);
%             y(2*M_.endo_nbr+order_var) = common+ghu*epsilon + abcOut2 + abcOut3;
%             resid(:,k) = irbc_small_so.sparse.dynamic_resid(y, x, M_.params, ys);
%         end
%         errors(:,t-1,n) = resid*int_weights;
%     end
% end
% toc;

