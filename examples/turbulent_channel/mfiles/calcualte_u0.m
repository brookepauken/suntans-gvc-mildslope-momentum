%% Calculate u0 

nu = 5e-4;
ustar = 0.02;
kappa = 0.4;
H = 10;

z0 = nu/(9*ustar);

u0 = ustar/kappa*(log(H/z0)+z0/H - 1)

Re_tau = ustar*H/nu

Nx = 400;

xspacing = 75/Nx * ustar/nu