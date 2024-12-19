%zeta(u)= u (u\le 1), 1 (1\le u \le 2) and 2 (u\ge 2)
function [riZ,riDZ,riD2Z]=irzetau(v,ep,zcase)
riZ=zeros(size(v));
% global zcase

%nargout is for the number of return values from [Zu,gradZ,D2Z]

if (nargout>1)
    riDZ=zeros(size(v));
    if (nargout>2)
        riD2Z=zeros(size(v));
    end
end
if (zcase==2)
%Zu:=Zeta(u)={u for u<=1, 1 for 1=<u<=2, u-1 for u>=2
riZ=v.*(v<=0)+(v/ep).*(v>0).*(v<=ep)+(v+1-ep).*(v>ep);
if (nargout>1)
    riDZ=(v<=0)+(1/ep).*(v>0).*(v<=ep)+(v>ep);
    if (nargout>2)
        %Diff Zeta''(u)=0 no-change in D2Z
    end

end
elseif (zcase==3)
    riZ=sqrt(v);
elseif (zcase==4)
    riZ=v;
end
