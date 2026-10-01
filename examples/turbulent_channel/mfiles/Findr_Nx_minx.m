function [rret, Nmax] = Findr_Nx_minx(length, Nx, minx, varargin)
if (nargin>3)
    max_r = varargin{1};
else
    max_r = 1.1;
end

syms r
rret = 0;
Nmax = Nx;
cnt = 0;
rolder = 0;
rold = 0;
nolder = Nx;
nold = Nx;

%dzbot = @(r) length*(r-1)/(r^(Nmax) - 1)-minx;

rval = fsolve(@(r)length*(r-1)/(r^(Nmax) - 1)-minx, max_r);
rret = rval;


% 
% 
% while rret > max_r | rret<1
%     if rret < 1
%         Nmax = Nmax-1;
%     elseif rret>max_r
%         Nmax = Nmax+1;
%     end
% 
%     if Nmax <= 1
%         if length/minx < 2
%             Nmax = 1;
%             rret = (length-minx)/minx;
%             break
%         else
%             Nmax = Nx+10;
%         end
% 
%     end
% 
% 
%     cnt= cnt+1;
%     assume(r,'real')
%     rval = vpasolve((length*(r-1)/(r^(Nmax) - 1))-minx, [0, Inf]);
% 
%     rret = eval(rval);
%     if isempty(rret)
%         rret = 0;
%     end
% 
%     if rolder == rret && rret ~= 0
%         rolds_new = [rolder, rold];
% 
%         %min_dis = min([max(rolds)-1.1, 1-min(rolds)]);
%         %rolds(min_dis == rolds-1.1)
%         if rret ~= min(nonzeros(rolds_new))
%             Nmax = nold;
%         end
%         rret=min(nonzeros(rolds_new));
% 
%         break
%     end
%     nolder = nold;
%     nold = Nmax;
%     rolder = rold;
%     rold = rret;
% 
%     if cnt>200;
%         rret = NaN;
%         break
% 
%     end
% end


end