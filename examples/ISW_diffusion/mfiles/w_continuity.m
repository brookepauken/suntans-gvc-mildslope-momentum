function [w_zedge, w_zcent, diffs] = w_continuity(w_djl, layer_thickness, u_edge, Nkmax, Nx)

%preallocate
w_zedge = zeros(Nkmax+1, Nx);
w_zcent = zeros(Nkmax, Nx);
diffs = zeros(Nkmax, Nx);
dx = 10000/Nx; %constant for now

%loop from bottom to top 

w_zedge(1, :) = 0; %bottom velocity is zero
for k = 2:(Nkmax+1)
    for i = 2:(Nx+1)
        dudx = (u_edge(k-1, i)-u_edge(k-1, i-1))/dx;
        dz = layer_thickness(k-1, i-1);
        w_zedge(k, i-1) = w_zedge(k-1, i-1) - dz*dudx;
        w_zcent(k-1, i-1) = .5*(w_zedge(k-1, i-1)+w_zedge(k, i-1));
        diffs(k-1, i-1) = w_djl(k-1, i-1) - w_zcent(k-1, i-1);
    end
end