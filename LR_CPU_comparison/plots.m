
% Select meshes
%meshes={'mesh_0016.mat';'mesh_0036.mat';'mesh_0064.mat';'mesh_0256.mat';'mesh_0400.mat';'mesh_1024.mat';'mesh_4096.mat'}
nbm=5;
F_Ndt=400;nbmeshes=4;
N=load(strcat('norms/T3_T_norms_bm',num2str(nbm),'RefM',num2str(F_Ndt),'nbmeshes',num2str(nbmeshes)));
h=[N.C_h;N.F_h]
L2zeta_error=N.Est_eznorm
H1zeta_error=N.GEst_eznorm
L1Xi_error=N.Est_eXi
 L2zeta=N.L2zeta;
 H1zeta=N.H1zeta;
 L1Xi=(N.L1Xi)

% er_Xi=N.er_Xi;
% er_znorm=N.er_znorm;
% Ger_znorm=N.er_znorm
%Ploting along spatial step h
% for imesh=1:nbmeshes-2 % convergence rate
%   ocL2zeta(imesh)=log(L2zeta_error(imesh)/L2zeta_error(imesh+1))/log(h(imesh)/h(imesh+1));
%   ocL1Xi(imesh)=log(L1Xi_error(imesh)/L1Xi_error(imesh+1))/log(h(imesh)/h(imesh+1));
%   ocH1zeta(imesh)=log(H1zeta_error(imesh)/H1zeta_error(imesh+1))/log(h(imesh)/h(imesh+1));
% end
% ploting along degree of freedom
% for imesh=1:nbmeshes-1 % convergence rate
%   ocL2p1_dof(imesh)=2*log(L2zeta_error(imesh)/L2zeta_error(imesh+1))/log(dof(imesh+1)/dof(imesh));
%   ocL1Xi_dof(imesh)=2*log(L1Xi_error(imesh)/L1Xi_error(imesh+1))/log(dof(imesh+1)/dof(imesh));
%   ocH1zeta_dof(imesh)=2*log(H1zeta_error(imesh)/H1zeta_error(imesh+1))/log(dof(imesh+1)/dof(imesh));
% end
% ocL2zeta
% ocL1Xi
% ocH1zeta

% p=loglog(h(1:size(h,1)-1),L1Xi_error,h(1:size(h,1)-1),L2zeta_error,h(1:size(h,1)-1),H1zeta_error,h(1:size(h,1)-1),h(1:size(h,1)-1));
 %axis([0.9*min(h(1:size(h,1)-1)) 1.1*max(h(1:size(h,1)-1)) 0.9*min([h(1:size(h,1)-1);L1Xi_error;L2zeta_error;H1zeta_error]) 1.1*max([h(1:size(h,1)-1);L1Xi_error;L2zeta_error;H1zeta_error])]);

%% Plot of the norms
  p=loglog(h,L1Xi,h,L2zeta,h,H1zeta);
  axis([0.8*min(h) 1.2*max(h) 0.8*min([L1Xi;L2zeta;H1zeta]) 1.2*max([L1Xi;L2zeta;H1zeta])]);
% error
%p=loglog(h(1:size(h,1)-1),L1Xi_error,h(1:size(h,1)-1),L2zeta_error,h(1:size(h,1)-1),H1zeta_error,h(1:size(h,1)-1),h(1:size(h,1)-1));
%axis([0.9*min(h(1:size(h,1)-1)) 1.1*max(h(1:size(h,1)-1)) 0.9*min([h(1:size(h,1)-1);L1Xi_error;L2zeta_error;H1zeta_error]) 1.1*max([h(1:size(h,1)-1);L1Xi_error;L2zeta_error;H1zeta_error])]);

%p=loglog(dof,L1Xi_error,dof,L2zeta_error,dof,H1zeta_error,dof,-dof);
%axis([0.9*min(dof) 1.1*max(dof) 0.9*min([-dof;L1Xi_error;L2zeta_error;H1zeta_error]) 1.1*max([-dof;L1Xi_error;L2zeta_error;H1zeta_error])]);

%%
% p=loglog(h(1:size(h,1)),Est_znorm_MT,h(1:size(h,1)-1),Est_znorm_HT);
% axis([0.9*min(h) 1.1*max(h) 0.9*min([Est_znorm_MT;Est_znorm_HT]) 1.1*max([Est_znorm_MT;Est_znorm_HT])]);
%%
%  p=loglog(h(1:size(h,1)-1),er_intxiu_MT,h(1:size(h,1)-2),er_intxiu_HT);
%  axis([0.9*min(h(1:size(h,1)-1)) 1.1*max(h(1:size(h,1)-1)) 0.6*min([er_intxiu_MT;er_intxiu_HT]) 1.1*max([er_intxiu_MT;er_intxiu_HT])]);
 %%
%p=loglog(h(1:size(h,1)),intxiu_MT,h(1:size(h,1)-1),intxiu_HT)
%axis([0.9*min(h) 1.1*max(h) 0.9*min([intxiu_MT;intxiu_HT]) 1.1*max([intxiu_MT;intxiu_HT])]);
grid on
% hold on
p(1).LineStyle=":";
p(2).LineStyle="--";
p(3).LineStyle="-."
p(1).Marker="o";
p(2).Marker="square";
p(3).Marker="*"
p(1).MarkerFaceColor="red";
p(2).MarkerFaceColor="blue";
p(3).MarkerFaceColor="black"
p(1).Color="red";
p(2).Color="blue";
p(3).Color="black"
p(1).LineWidth=1.5;
p(2).LineWidth=1.5;
p(3).LineWidth=1.5
%%
  lgdEZ=legend({'$E^{m}_{\Pi_{\mathcal{D}}\Xi}$','$E^{m}_{\Pi_{\mathcal{D}}\zeta}$','$E^{m}_{\nabla_{\mathcal{D}}\zeta}$','Ref_slope one'},'Location','southeast','Interpreter','latex');
 title(lgdEZ);

%  lgdEZ=legend({'$\Pi_{\mathcal{D}}\Xi$','$\Pi_{\mathcal{D}}\zeta$','$\nabla_{\mathcal{D}}\zeta$','Ref_slope one'},'Location','southeast','Interpreter','latex');
%  title(lgdEZ);
 
 %lgdEZ=legend({'L1Xi_error','L2zeta_error','H1zeta_error','Ref_slope'},'Location','northeast');%,'Interpreter','latex');
 %title(lgdEZ,"MLP1: various norms wrt dof" ,'Interpreter','latex');
 
 %%
% lgdCZ=legend({'MLP1-T','HMM-T'},'Location','northwest');
% title(lgdCZ,"Convergence of $\|\zeta(u)\|_{L^{2}(\Omega\times (0,T)\times \Theta)}$",'Interpreter','latex')
%%
% lgdEx=legend({'MLP1-T','HMM-T'},'Location','southeast','Interpreter','latex');
%  title(lgdEx,"Comparison plot of $E^{m}_{\Xi}$" ,'Interpreter','latex');
 %%
% lgdCx=legend({'MLP1-T','HMM-T'},'Location','northwest');
%title(lgdCx,"Convergence of $\|\Xi(u_{T})\|_{L^{1}(\Omega\times \Theta)}$",'Interpreter','latex')