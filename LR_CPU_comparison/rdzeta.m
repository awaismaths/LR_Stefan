%zeta(u)= u (u\le 1), 1 (1\le u \le 2) and 2 (u\ge 2)
function rd=rdzeta(u,ep,zcase)
rd=zeros(size(u));
% global zcase
if (zcase==2)
rd=(u<0)+(ep).*(u>0).*(u<=1)+(u>1);
end
end

