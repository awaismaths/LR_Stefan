function [ITER,Res]=timestepping_RLSSP(tol,ep,Lc,zcase,ucase,bmm,dW,ncell,nvert,area,cg,cell_v,X,Mass,L,A,G,B_indices,I_indices,imesh,Ndt,dt,bdry_vert,mesh,h)
%global ucase zcase
itermax=5000;
usol_idt=zeros(nvert,Ndt+1);
usol_idt(:,1)=X;
%fid = fopen('results.txt','w');
ITER=0;
%num_updates=0;
Res=0;
zb=zeros(nvert,1);
for idt=1:Ndt;
    %str = sprintf('idt=%d / %d\n',idt,Ndt);
    % forkprint(fid,str);
    Xprev = X;
    %Source: Compute using (u_t)-laplace(zeta) and the exact solution
    b=assemble_source_stc(idt,cell_v,ncell,nvert,area,dW,Xprev,zcase);
    % Following is the source function for the deterministic SP
    % constructed from the exact solution.
    %b=assemble_source(idt*dt,cell_v,ncell,nvert,area,cg,zcase,ucase);
    %Dirichlet non-homogeneous BC
    %zb(B_indices)=-1*ones(size(B_indices,1),1);
    zb(B_indices)=zetau((exact_solution(idt*dt,bdry_vert(B_indices,:),ucase))',zcase);%Zb(:,idt);%
    rhs = Mass*Xprev + dt*b;
    iter = 0;
    res = 1;
    while (iter < itermax && res > tol)
        B=-dt*A*(rzetau(Xprev,ep,zcase)-Lc*Xprev)+rhs;
        X_b=zeros(nvert,1);%Xprev;
        X_b(B_indices)=irzetau(zb(B_indices),ep,zcase);%(1/Lc)*(zb(B_indices)-zetau(Xprev(B_indices),zcase))+Xprev(B_indices);
        F=B-G*X_b;
        % Using Cholesky decomposition
        Xsol=L\(L'\F(I_indices));
        X=X_b;
        X(I_indices)=Xsol;
        iter = iter+1;
        % residual by increments
        %res = norm(X-Xprev,Inf) / norm(Xprev,Inf);
        Xprev=X;
        % residual of non-linear system
        residue=Mass*X+dt*A*rzetau(X,ep,zcase)-rhs;
        res = norm(residue(I_indices), 2);
    end; % end nonlinear iterations
    if (iter==itermax)
        res
        iter
        error('no convergence')
    end;
    usol = X;
     usol(B_indices)=rzetau(X(B_indices),ep,zcase);
     usol_idt(:,idt+1)=usol;
    ITER=ITER+iter;
    Res=Res+abs(res);
    % ex_sol=exact_solution(idt*dt,vertex)';
    %Res = norm(ex_sol(B_indices)-X(B_indices),Inf) / norm(ex_sol(B_indices),Inf)

    % Write the solution and grid vtk files, to be plotted by "paraview"
    %         write_solution_vtk(usol,strcat('VTKout/solution',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);

    % Exact solution
    %     ex_sol=exact_solution(dt*idt,vertex,ucase)';
    %         write_solution_vtk(ex_sol,strcat('VTKout/ex_sol',num2str(idt)),ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
    %str = sprintf('Solution computed, iter=%d, res=%4.2e, max sol=%f, max ex_sol=%f\n', iter, res, max(usol(I_indices)), max(ex_sol(I_indices)));
    %forkprint(fid,str);
    %z_error=norm(zetau(usol(B_indices))-zetau(ex_sol(B_indices)),"inf")
end; % end time stepping
ITER=(ITER/Ndt);
Res=Res/Ndt;
%num_updates=ceil(num_updates/Ndt);
%str = sprintf('Solution computed for bm=%d,num_time_steps=%d, Avg_iter=%d Avg_res=%4.2e\n', bmm,Ndt,ITER,Res);
%forkprint(fid,str);
%Saving solutions for each Brownian motion
save(strcat('solutions/RBM',num2str(bmm),'mesh',mesh(1:8),'tcuz',num2str(ucase),num2str(zcase)),'usol_idt','dt','Ndt','mesh','h','-v7.3');
end