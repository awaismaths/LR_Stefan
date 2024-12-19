% Assemble the global diffusion matrix
function F = assemble_source_stc(idt,cell_v,ncell,nvert,area,dW,Xprev,zcase)
F=zeros(nvert,1);
%global zcase
    for i=1:ncell
        if (zcase==2)
        F(cell_v{i}(1:3))=F(cell_v{i}(1:3))+500*(area(i)/3)*dW(idt)*sqrt(Xi(Xprev(cell_v{i}(1:3)),zcase));
     elseif (zcase==1)
         F(cell_v{i}(1:3))=F(cell_v{i}(1:3))+(area(i)/3)*dW(idt)*sqrt(Xi(Xprev(cell_v{i}(1:3))+ones(3,1),zcase)-Xprev(cell_v{i}(1:3)));
     end
        %F(cell_v{i}(1:3))=F(cell_v{i}(1:3))+(area(i)/3)*dW(idt)*sqrt(Xi(Xprev(cell_v{i}(1:3)),zcase));
    
    end
end