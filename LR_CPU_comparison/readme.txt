%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %%%% CPU TIME COMPARISON %%%%
% Mass-lumped P1 codes with different solvers
% for u_t-Delta zeta(u)=0 with Dirichlet BC
%
%    Author: Muhammad Awais Khan
%    Date: 19/12/2024
% % %
% This code is provided as a companion of the article
%   "Efficient iterative linearised solvers for numerical approximations of stochastic Stefan problems",
%
%
% Usage and modification of this code is permitted, but any scientific publication resulting from
% part of this code should cite the aforementioned article.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 This code is used to compare the CPU time for Newton, Linearised and Regularised solvers for deterministic and stochastic case
Run the files in the following order: 

For stochastic case, use "main_Smin.m" and change the parameters C_e for ep (regularisation constant) and tolerance parameter c_tol. These parameters can be selected from {1,10,100,1000}
For deterministic case, use "main_dmin.m" and use the same above instruction. Look for the "dat" files to see the CPU time results after if code runs successfully.

You need to run "gen_bm.m" once when you first run the code to generate and save a sample of Brownian motions. Select the number of Brownian motions and appropriate time steps for the selected meshes.
