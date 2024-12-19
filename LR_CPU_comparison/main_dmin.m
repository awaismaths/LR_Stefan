%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPU Time Comparison
% Mass-lumped P1 codes with different solvers (Newton, Linearised,
% Regularised)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long
T=1;
dt_initial=[0.1;0.01];%0.001;0.0001]
%Lc_M=[0.64 0.60 0.75 0.87;1 0.58 0.61 0.75;0.53 0.53 0.56 0.65;0.52 0.51 0.52 0.57];
% Lc_M=[0.64	0.59	0.75	0.87 1;
% 1	0.58	0.61	0.75 0.89;
% 0.53	0.53	0.56	0.65 0.76;
% 0.52	0.51	0.52	0.57 0.67;
% 0.51	0.51	0.51	0.53 0.58];

nts=size(dt_initial,1);
meshes={'mesh1_1.mat';'mesh1_2.mat'};%'mesh1_3.mat';'mesh1_4.mat';'mesh1_5.mat';'mesh1_6.mat'};%'mesh1_7.mat'};
nbmeshes=size(meshes,1);
% % Write data file
fid = fopen('dmin.dat','w');
fprintf(fid, 'sr h Ndt NI LI RI NR LR RR NC NCC LC LCC RC RCC RNL RNR\n');
t=tiledlayout(nbmeshes,2);
xlabel(t,'dt');
ylabel(t,'CPU Time');
axis='padded';
t.Padding = 'compact';
t.TileSpacing = 'compact';
MTime_N=zeros(nbmeshes,nts);
MTime_L=zeros(nbmeshes,nts);
MTime_R=zeros(nbmeshes,nts);

% Calculate cumulative frequency for time
MTN = 0;
MTL = 0;
MTR = 0;
sr = 0;

% Loop over each mesh
for imesh = 1:nbmeshes
    load(strcat('../HHO-Lapl-OM-master/matlab_meshes/',meshes{imesh}));
     h=max(abs(diam));
    for j=1:nts;        
         h;
         dt_initial(j);
        %ep=max(dt_initial(j)^4,h^4)%min(h,dt_initial(j))^2
        %%%%%%%%%%%%%
        %%%%% Set your tolerance and regularisation parameters here!!!
        C_e=1;
        C_tol=10;
        tol=(min(dt_initial(j),h))/C_tol;%max((dt_initial(j))^2,h^2)/1e4
        ep=min(dt_initial(j),h)/C_e;
         for k=1:3;
            [ITER_N(k),TN(k),Res_N(k)]=sp(tol,nvert,cell_n,diam,ncell,vertex,area,cell_v,dt_initial(j),meshes{imesh});
            [ITER_L(k),TL(k),Res_L(k)]=lsp(tol,0.63,nvert,cell_n,diam,ncell,vertex,area,cell_v,dt_initial(j),meshes{imesh});
            [ITER_R(k),TR(k),Res_R(k)]=rlsp(ep,tol,0.63,nvert,cell_n,diam,ncell,vertex,area,cell_v,dt_initial(j),meshes{imesh});
            disp([imesh,j,k])
            %Lc_M(imesh,j)
        end
        sr=sr+1;
        MTime_N(imesh,j)=min(TN);
        MTime_L(imesh,j)=min(TL);
        MTime_R(imesh,j)=min(TR);
        Ratio_LN=(MTime_N(imesh,j))/MTime_L(imesh,j);%((MTime_L(imesh,j)-MTime_N(imesh,j))/MTime_L(imesh,j))*100;
        Ratio_RN=(MTime_N(imesh,j))/MTime_R(imesh,j);%((MTime_R(imesh,j)-MTime_N(imesh,j))/MTime_R(imesh,j))*100;
        MTN = MTN + MTime_N(imesh, j);
        MTL = MTL + MTime_L(imesh, j);
        MTR = MTR + MTime_R(imesh, j);
       
        fprintf(fid,'%d %f %d %.2f %.2f %.2f %.2e %.2e %.2e %.4f %.4f %.4f %.4f %.4f %.4f %.2f %.2f\n',...
            sr,h,T/dt_initial(j), min(ITER_N),min(ITER_L),min(ITER_R),mean(Res_N),...
            mean(Res_L),mean(Res_R),min(TN), MTN, min(TL),MTL, min(TR),MTR, Ratio_LN, Ratio_RN);
        %ATime_N(imesh,j),ATime_L(imesh,j),percentageA(imesh,j)
    end
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
exportgraphics(t, 'CompPlotD.jpg');

% Close the data file
fclose(fid);