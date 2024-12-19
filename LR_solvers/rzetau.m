%zeta(u)= u (u\le 1), 1 (1\le u \le 2) and 2 (u\ge 2)
function [rZu,rDZ,rD2Z]=rzetau(u,ep,zcase)
rZu=zeros(size(u));
% global zcase

%nargout is for the number of return values from [Zu,gradZ,D2Z]

if (nargout>1)
    rDZ=zeros(size(u));
    if (nargout>2)
        rD2Z=zeros(size(u));
    end
end
if (zcase==2)
    %Zu:=Zeta(u)={u for u<=1, 1 for 1=<u<=2, u-1 for u>=2
    rZu=u.*(u<=0)+(ep*u).*(u>0).*(u<=1)+(u-1+ep).*(u>1);
    if (nargout>1)
        rDZ=(u<0)+(ep).*(u>0).*(u<=1)+(u>1);
        if (nargout>2)
            %Diff Zeta''(u)=0 no-change in D2Z
        end

    end
elseif (zcase==3)
    rZu=u.^(2);
elseif (zcase==4)
    rZu=u;
end
end

