%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Run the files in the following order: 

* Step 1: Run "MAIN_MLP1_*****.m" (In place of *****, for the stochastic case: choose NSSP for Newton, LSSP for linearised, RLSSP for regularised solver.
	  While in the deterministic case: choose NSP for Newton, LSP for linearised, RLSP for regularised solver.
* Step 2 [see also Note 1 below]: Run "compare_time.m", "gen_bm" and "compare_meshes.m" to prepare the comparison between norms of solutions.
* Step 3: Run "Compute_norms_*****.m" to compute the L2 error of: \zeta(u) and \nabla\zeta(u) and L1 error of \Xi(u) (see step 1 to fill *****).

Note 1: Step 2 only has to be run once for each mesh family; the files created in the Time-comparisons, BM and Mesh-comparisons folders
during this step can be re-used for any solution computed on the same mesh families. You also need to run "gen_bm.m" once when you first run the code to generate and save a sample of Brownian motions. Select the number of Brownian motions and appropriate time steps for the selected meshes.
Note 2: The averaged solutions can be generated by "writeVtkforBM.m" and viewed in Paraview.
Note 3: Deterministic case does not require step 2 and step 3.
