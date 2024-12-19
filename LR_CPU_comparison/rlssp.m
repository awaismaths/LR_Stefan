function [ITER_RLSSP,TRLSSP,Res_RLSSP]=main_comp_rlssp(tol,ep,Lc,nbm,nvert,cell_n,diam,ncell,vertex,area,cell_v,dt,mesh)
format long;
imesh=1;
meshes{imesh}=mesh;
format long;
%global ucase zcase;
ucase=2;
zcase=2;
% Final times
T=1;
%%
% Sequence of meshes over which we want to run the scheme
%%
% The meshes are available at https://github.com/jdroniou/HHO-Lapl-OM
%meshes={'mesh1_1.mat';'mesh1_2.mat';'mesh1_3.mat'};%'mesh1_4.mat'};%'mesh1_5.mat'};%'mesh1_06.mat'};
nbmeshes=size(meshes,1);
% Lp1_error=zeros(nbmeshes,1);
% L2zeta_error=zeros(nbmeshes,1);
h=zeros(nbmeshes,1);
Ndt=zeros(nbmeshes,1);
% To see the results printed in file
fid = fopen('results.txt','w');
%%%fclose(fid);
Ndt(1) = ceil(T/dt);
for imesh=1:nbmeshes
    % Load mesh
    %loadmesh=strcat('load ../../HHO-Lapl-OM-master/matlab_meshes/',meshes{imesh});
    %str = sprintf('%s\n',loadmesh);
    % forkprint(fid,str);
    %eval(loadmesh);
    % Compute real centers of gravity, mesh diameter and area of dual mesh
    cg=gravity_centers(ncell,cell_v,vertex,area);
    h(imesh)=max(abs(diam));%
    %% Finding boundary vertices
    fbc=zeros(size(vertex,1),1);
    bdry_vert=zeros(size(vertex,1),2);
    for i=1:ncell
        I=find(cell_n{i}==0);
        if (size(I,2)>0)
            bdry_vert_indices = [cell_v{i}(I) cell_v{i}(I+1)];
            bdry_vert(bdry_vert_indices,:)  = vertex(bdry_vert_indices,:);
            fbc(bdry_vert_indices)=1;
        end
    end
    I_indices=find(~fbc);
    B_indices=find(fbc);
    %nbvert=size(B_indices,1);
    dualarea=compute_dualarea(area,ncell,nvert,cell_v,B_indices);
    % str = sprintf('mesh= %s, h= %4.2e, time step= %4.2e \n',meshes{imesh},h(imesh),T/Ndt(imesh));
    %forkprint(fid,str);

    %% Initialise RHS and unknown
    % Initial condition and exact solution
    ex_sol=exact_solution(0,vertex,ucase)'; % Exact solution at t0=0

    %     write_solution_vtk(ex_sol,strcat('VTKout/solution0'),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    %     write_solution_vtk(ex_sol,'VTKout/ex_sol0',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    % Unknows at the boundary are "zeta(u)" while "u" at the interior points
    X=ex_sol;
    % X(B_indices)=zetau(ex_sol(B_indices),zcase);
    %% Loading brownian motion
    %%%%Constructing Brownian motion (Finacial toolbox is  required)
%     W=zeros(Ndt(imesh),nbm);
%     for bmm=1:nbm
%         [Path,Time,dW]=simulate(bm(0,1),Ndt(imesh));
%         W(:,bmm)=dW;
%     end
    W=zeros(Ndt(imesh),nbm);
    load(strcat('BM/BM_mesh',num2str(meshes{nbmeshes}(1:7)),'Ndt',num2str(Ndt(imesh))),'WF');
    WF=WF(:,1:nbm);
      if imesh<nbmeshes
        for i=1:Ndt(imesh)
            db=floor(Ndt(nbmeshes)/Ndt(imesh))
            W(i,:)=sum(WF([(db*(i-1)+1):(i*db)],:),1);

        end
        else
        W=WF;
     end
    clear WF
    %% Time starts here !!!
    tstartR=tic;
    %% ASSEMBLE MATRIX of Laplacian
    [A,b]=assemble_diffusion_system(cell_v,ncell,nvert,vertex);
    Mass = spdiags(dualarea,0,nvert,nvert);
    G=Mass+dt*Lc*A;
    L=chol(G(I_indices,I_indices));
    %     Zb=zeros(nbvert,Ndt(imesh));
    %     for idt=1:Ndt(imesh)
    %          %dualarea(BC_indices).*(exact_solution(idt*dt,bdry_vert(BC_indices,:),ucase)) ...
    %           %  + forcediff*dt.*zetau((exact_solution(idt*dt,bdry_vert(BC_indices,:),ucase)),zcase);
    %         Zb(:,idt)=zetau((exact_solution(idt*dt,bdry_vert(B_indices,:),ucase))',zcase);
    %     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for bmm=1:nbm
        [iter,res]=timestepping_RLSSP(tol,ep,Lc,zcase,ucase,bmm,W(:,bmm),ncell,nvert,area,cg,cell_v,X,Mass,L,A,G,B_indices,I_indices,imesh,Ndt(imesh),dt,bdry_vert,meshes{imesh},h);
        ITER_RLSSP(bmm)=iter;
        Res_RLSSP(bmm)=res;
%         % Following error only works in the case of deterministic problem. Also, uncomment the correct source in the "timestepping".
%         Sol=load(strcat('solutions/RLSSP_BM',num2str(bmm),'mesh',meshes{imesh}(5:7),'tcuz',num2str(ucase),num2str(zcase)));
%         Sol=Sol.usol_idt;
%         % Exact solution
%         ex_sol=exact_solution(dt*Ndt(imesh),vertex,ucase)';
%         usol=Sol(:,Ndt(imesh)+1);
%         
%         % compute error at the final times
%         [Lp1_error(imesh) L2zeta_error(imesh) H1zeta_error(imesh)] = R_compute_errors(usol,ex_sol,dualarea,A,I_indices,zcase,ucase,ep);
%         str = sprintf('Mesh %i. Errors: L^(2)=%4.2e, L^2 on zeta(u):%4.2e, H1 on zeta(u):%4.2e\n',imesh,Lp1_error(imesh),L2zeta_error(imesh),H1zeta_error(imesh));
%         forkprint(fid,str);
    end
end; % end meshes
Time=toc(tstartR);
ITER_RLSSP=mean(ITER_RLSSP);
Res_RLSSP=mean(Res_RLSSP);
TRLSSP=Time;
%disp(['Elapsed time is ' num2str(toc) 'seconds' ])

end






