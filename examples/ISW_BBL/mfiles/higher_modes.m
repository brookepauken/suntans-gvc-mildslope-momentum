%
% HIGHER_MODES
%
% [w,c]=higher_modes(Z,N,OMEGA,NMODE,TOPBC) computes the first
% NMODE modes of solution of the eigenvalue problem
%
% w_zz + (N^2-omega^2)/c^2 w = 0
%
% subject to the no-flux boundary condition at the bottom w=0, and
% at the free surface:
%   1) TOPBC = 'rigid lid': The first mode is the first baroclinic mode
%   2) TOPBC = 'free surface': The first mode is the first barotropic mode
%
% and where:
%   N(z) is the buoyancy frequency
%   omega is the wave frequency (note that omega<N)
%   c is the wave speed. When written as
%
% 1/(N^2-omega^2) w_zz + 1/c^2 w = 0
%
% we have the eigenvalue problem where w is the eigenfunction corresponding
% to the eigenvalue lambda = -1/c^2.
%
% For more more information, refer to the appendix in Fringer and Street (2003)
% doi:10.1017/S0022112003006189
%

%
% Oliver Fringer
% Stanford University
%
function [w,c] = higher_modes(z,N,omega,nmode,topbc)

     g = 9.81;
     if(diff(z)<0)
       error('z must be monotonically increasing');
     end

     Nk = length(N);
     dz = zeros(Nk,1);
     dziph = zeros(Nk+1,1);

     dziph(2:end-1) = z(2:end)-z(1:end-1);
     dziph(1) = 2*dziph(2)-dziph(3);
     dziph(end) = 2*dziph(end-1)-dziph(end-2);

     % a(k)*w(k-1)+b(k)*w(k)+c(k)*w(k+1) + f/c^2*w(k) = 0     
     a = zeros(Nk,1);
     b = zeros(Nk,1);
     c = zeros(Nk,1);
    
     dz = 0.5*(dziph(1:end-1)+dziph(2:end));
     a = 1./(dziph(1:end-1).*dz);
     b = -(1./dziph(1:end-1)+1./dziph(2:end))./dz;
     c = 1./(dziph(2:end).*dz);

     f = (N.^2-omega^2);
     a = -a./f;     
     b = -b./f;
     c = -c./f;

     % Bottom no-flux w=0.
     b(1) = b(1)-a(1);

     if(strcmp(topbc,'rigid lid'))
       % Top rigid lid: w=0 
       b(Nk) = b(Nk)-c(Nk);
     elseif(strcmp(topbc,'free surface'))
       % Free-surface: dw/dz + g/c^2 w = 0
       a(Nk) = -1/(g*dziph(end-1));
       b(Nk) = 1/(g*dziph(end-1));
     else
       error('topbc input must be either ''rigid lid'' or ''free surface''');
     end

     A = diag(a(2:end),-1)+diag(b,0)+diag(c(1:end-1),1);
     [X,L] = eig(A,'nobalance');

     [lams,is] = sort(diag(L));
     X = X(:,is);

     w = X(:,[1:nmode]);
     c = sqrt(1./lams([1:nmode]));

