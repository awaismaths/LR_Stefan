
%dzeta(u)= 1 (u< 1), 0 (1\le u \le 2) and 1 (u> 2)
function dzu=dzeta(u,zcase)
dzu=zeros(size(u));
% global zcase
if (zcase==1)
        dzu=(u<1)+(u>2);
    elseif (zcase==2)
       dzu=(u<0)+(u>1);
    elseif (zcase==3)
    dzu=2*u;
    elseif (zcase==4)
    dzu=1;
end
end
