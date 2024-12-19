%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPU Time Comparison
% Mass-lumped P1 codes with different solvers (Newton, Linearised,
% Regularised)
% for u_t-Delta zeta(u)=f(zeta(u))*dW with Dirichlet BC
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear all variables and set format to long for better precision
clear all
format long

% Parameters
T = 1; % Total time
nbm = 10; % Number of Brownian motions
Lc = 0.63; % Lipschitz constant
dt_initial = [0.1; 0.01];% 0.001; 0.0001]; % Initial time steps

divd = [1e1];% 1e2]; % Division factors for epsilon or tol
ndivd = size(divd, 1); % Number of division factors
nts = size(dt_initial, 1); % Number of time steps
meshes = {'mesh1_1.mat'};% 'mesh1_2.mat'; 'mesh1_3.mat'};% 'mesh1_4.mat'; 'mesh1_5.mat'};% 'mesh1_6.mat'}; % Mesh files
nbmeshes = size(meshes, 1); % Number of meshes

% Open data file for writing
fid = fopen('Smin.dat', 'w');
fprintf(fid, 'sr h Ndt NI LI RI NR LR RR NC NCC LC LCC RC RCC RNL RNR\n');

% Prepare for tiled layout for plots
t = tiledlayout(nbmeshes, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
xlabel(t, 'dt');
ylabel(t, 'CPU Time');

% Calculate cumulative frequency for time
MTN = 0;
MTL = 0;
MTR = 0;
sr = 0;
% Loop over each mesh
for imesh = 1:nbmeshes
    % Load mesh data
    load(strcat('../HHO-Lapl-OM-master/matlab_meshes/', meshes{imesh}));
    h = max(abs(diam));
    Ndt(imesh)=2*round(0.5*T/h^2);
    dt=T/Ndt(imesh);

    % Initialize arrays to store minimum times
    MTime_N = zeros(nbmeshes, nts);
    MTime_L = zeros(nbmeshes, nts);
    MTime_R = zeros(nbmeshes, nts);

    % Loop over each time step
    for j = 1:nts
        %dt = dt_initial(j);
        %%%%%%%%%%%%%
        %%%%% Set your tolerance and regularisation parameters here!!!
        C_e=1;
        C_tol=10;
        ep = min(dt, h) / C_e;
        tol = min(dt, h) / C_tol;
        
        % Initialize arrays to store results for each iteration
        rep=1;
        TN = zeros(1, rep);
        TL = zeros(1, rep);
        TR = zeros(1, rep);
        ITER_N = zeros(1, rep);
        ITER_L = zeros(1, rep);
        ITER_R = zeros(1, rep);
        Res_N = zeros(1, rep);
        Res_L = zeros(1, rep);
        Res_R = zeros(1, rep);

        % Perform rep iterations for each configuration
        for k = 1:rep
            [ITER_N(k), TN(k), Res_N(k)] = ssp(tol, nbm, nvert, cell_n, diam, ncell, vertex, area, cell_v, dt, meshes{imesh});
            [ITER_L(k), TL(k), Res_L(k)] = lssp(tol, 0.63, nbm, nvert, cell_n, diam, ncell, vertex, area, cell_v, dt, meshes{imesh});
            [ITER_R(k), TR(k), Res_R(k)] = rlssp(tol, ep, 0.63, nbm, nvert, cell_n, diam, ncell, vertex, area, cell_v, dt, meshes{imesh});
            disp([imesh, j, k])
        end

        % Store minimum times and calculate ratios
        MTime_N(imesh, j) = min(TN);
        MTime_L(imesh, j) = min(TL);
        MTime_R(imesh, j) = min(TR);
        sr = sr + 1;

        MTN = MTN + MTime_N(imesh, j);
        MTL = MTL + MTime_L(imesh, j);
        MTR = MTR + MTime_R(imesh, j);
        Ratio_LN = min(TN)/min(TL);
        Ratio_RN= min(TN) / min(TR);

        % Write results to file
        fprintf(fid, '%d %f %d %.2f %.2f %.2f %.2e %.2e %.2e %.4f %.4f %.4f %.4f %.4f %.4f %.2f %.2f\n',...
            sr,h, T / dt, min(ITER_N), min(ITER_L), min(ITER_R), mean(Res_N), mean(Res_L),...
            mean(Res_R), min(TN), MTN, min(TL),MTL, min(TR),MTR, Ratio_LN, Ratio_RN);
    end

    % Plot results for the current mesh
    nexttile;
    loglog(dt_initial, MTime_N(imesh, :), '-o', 'DisplayName', 'Newton', 'LineWidth', 1.5);
    hold on;
    loglog(dt_initial, MTime_L(imesh, :), '-s', 'DisplayName', 'Linearised', 'LineWidth', 1.5);
    hold on;
    loglog(dt_initial, MTime_R(imesh, :), '-d', 'DisplayName', 'Regularised', 'LineWidth', 1.5);
    hold off;

    title(['Mesh: ', meshes{imesh}(1:8), ' with minimum time']);
    xlabel('Time Step (dt)');
    ylabel('CPU Time');
    legend('show', 'Location', 'best');
    grid on;
  
end

  % Add legend to the plot
lgd = legend('Newton', 'Linearised', 'Regularised');
lgd.Layout.Tile = 'south';

% Save the plot
exportgraphics(t, 'CompPlotS.jpg');

% Close the data file
fclose(fid);