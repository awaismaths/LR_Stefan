clear all;

format long;
tic
% Parameters for nonlinearity and barenblatt solution
global t0;
global CB;
%global ucase zcase;
ucase=2;
zcase=2;

forcediff=1;

t0=0.1;
CB=0.005;
%t0=.5;
%CB=0.1;

% Final times
T=1;
dt_initial=0.25;%0.001;

% nonlinear iterations
tol=1e-12;
itermax=1000;
% Relaxation parameter (1 for no relaxation)
relax = 1;

%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat';'mesh1_4.mat';'mesh1_5.mat';'mesh1_6.mat'};%'mesh1_7.mat'};

nbmeshes=size(meshes,1);
Lp1_error=zeros(nbmeshes,1);
L2zeta_error=zeros(nbmeshes,1);
h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);

% To see the results printed in file
fid = fopen('results.txt','w');
str = sprintf('t0=%f\n',t0);%sprintf('m=%f, t0=%f\n',m,t0);
forkprint(fid,str);
%%%fclose(fid);
Ndt(1) = ceil(T/dt_initial);

for imesh=1:nbmeshes
    % Load mesh
    loadmesh=strcat('load ../HHO-Lapl-OM-master/matlab_meshes/',meshes{imesh});
    str = sprintf('%s\n',loadmesh);
    forkprint(fid,str);
    eval(loadmesh);
    % Compute real centers of gravity, mesh diameter and area of dual mesh
    cg=gravity_centers(ncell,cell_v,vertex,area);
    h(imesh)=max(abs(diam));%
   numcells(imesh)=ncell;
   numedges(imesh)=nedge;
   numvert(imesh)=nvert;
end; % end meshes
h
numcells
numedges
numvert