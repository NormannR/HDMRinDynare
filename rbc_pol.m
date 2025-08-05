addpath /home/normann/Dropbox/dynare/matlab;

dynare rbc_so;
%%
dr = oo_.dr;
ys = dr.ys;
k2 = M_.nstatic+(1:M_.nspred);
order_var = dr.order_var;
order_var_k2 = order_var(k2);
ys_order_var = ys(order_var);
ys_order_var_k2 = ys(order_var_k2);
%%
ind_c = find(strcmp(M_.endo_names, 'c'));
row_c = dr.inv_order_var(ind_c);
ind_z = find(strcmp(M_.endo_names, 'z'));
col_z = find(ind_z == dr.state_var);
ind_k = find(strcmp(M_.endo_names, 'k'));
row_k = dr.inv_order_var(ind_k);
col_k = find(ind_k == dr.state_var);
ind_rk = find(strcmp(M_.endo_names, 'rk'));
row_rk = dr.inv_order_var(ind_rk);
%%
rho = M_.params(strcmp(M_.param_names, "rho"));
sig = M_.params(strcmp(M_.param_names, "sig"));
%%
c = ys(row_c);
c_k = dr.ghx(row_c,col_k);
c_z = dr.ghx(row_c,col_z)/rho;
% Check
-dr.ghx(row_c,col_z)*sig/rho+dr.ghu(row_c,1)
%%
c_kk = 0.5*dr.ghxx(row_c,col_k);
c_kz = dr.ghxx(row_c,col_z)/rho;
c_zz = 0.5*dr.ghxx(row_c,2*col_z)/rho^2;
c_0 = 0.5*dr.ghs2(row_c)+0.5*dr.ghuu(row_c)+0.5*dr.ghxx(row_c,2*ind_z)*sig^2/rho^2;
% Check
dr.ghxu(row_c,ind_z)-dr.ghxx(row_c,2*col_z)*sig/rho
dr.ghxu(row_c,ind_k)-dr.ghxx(row_c,col_z)*sig/rho
%%
rk = ys(row_rk);
rk_k = dr.ghx(row_rk,col_k);
rk_z = dr.ghx(row_rk,col_z)/rho;
% Check
-dr.ghx(row_rk,col_z)*sig/rho+dr.ghu(row_rk,1)
%%
rk_kk = 0.5*dr.ghxx(row_rk,col_k);
rk_kz = dr.ghxx(row_rk,col_z)/rho;
rk_zz = 0.5*dr.ghxx(row_rk,2*col_z)/rho^2;
rk_0 = 0.5*dr.ghs2(row_rk)+0.5*dr.ghuu(row_rk)+0.5*dr.ghxx(row_rk,2*ind_z)*sig^2/rho^2;
% Check
dr.ghxu(row_rk,ind_z)-dr.ghxx(row_rk,2*col_z)*sig/rho
dr.ghxu(row_rk,ind_k)-dr.ghxx(row_rk,col_z)*sig/rho
%%
A_k = dr.ghx(:,col_k);
A_z = dr.ghx(:,col_z)/rho;
g_1 = [A_k A_z];
% Check
-dr.ghx(:,col_z)*sig/rho+dr.ghu(:,1)
%%
A_kk = 0.5*dr.ghxx(:,col_k);
A_kz = 0.5*dr.ghxx(:,col_z)/rho;
A_zz = 0.5*dr.ghxx(:,2*col_z)/rho^2;
g_2 = [A_kk A_kz A_kz A_zz];
g_0 = ys_order_var+0.5*dr.ghs2+0.5*dr.ghuu+0.5*dr.ghxx(:,2*ind_z)*sig^2/rho^2;
% Check
dr.ghxu(:,ind_z)-dr.ghxx(:,2*col_z)*sig/rho
dr.ghxu(:,ind_k)-dr.ghxx(:,col_z)*sig/rho
%%
save rbc_pol.mat row_c row_k row_rk order_var order_var_k2 ys_order_var ys_order_var_k2 g_0 g_1 g_2;