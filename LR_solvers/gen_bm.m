%% Brownian motions for the finest mesh to be used for all other coarse meshes
%Reminder: we are using dt=h^2;
%Constructing Brownian motion (Finacial toolbox is  required)
clear all
T=1;
nbm=500;
meshes={'mesh1_11.mat';'mesh1_12.mat';'mesh1_13.mat';'mesh1_14.mat';'mesh1_15.mat';'mesh1_16.mat'};
nbmeshes=size(meshes,1);
for imesh=1:nbmeshes
fine_mesh=matfile(strcat('../matlab_meshes/',meshes{imesh}));
h(imesh)=max(abs(fine_mesh.diam));
Ndt(imesh)=2*round(0.5*T/h(imesh)^2);
dt=T/Ndt(imesh);
WF=zeros(Ndt(imesh),nbm);
     for bmm=1:nbm
         [Path,Time,dW]=simulate(bm(0,1),Ndt(imesh));
         WF(:,bmm)=sqrt(dt)*dW;
     end
     save(strcat('BM/BM_mesh',num2str(meshes{imesh}(1:8))),'WF');
     clear WF;
end