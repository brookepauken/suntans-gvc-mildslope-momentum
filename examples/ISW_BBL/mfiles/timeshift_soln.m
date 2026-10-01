function u_interp = timeshift_soln(c, dt, L, numshift, X, Z, u, Z_newgrid)
%shift quantity u back c*dt numshift times on periodic grid of length L

dis_shift = numshift*c*dt;
[r, ~] = size(X);
dx = X(1, 2)-X(1, 1);

X_shift = X - dis_shift; %shift back how far wave propogated in numshift time steps
[~, inds_outofbox] = find(X_shift<=0);
inds_outofbox=unique(inds_outofbox);
if isempty(inds_outofbox)
    X_shift = [X_shift, ones(r, 1)*(L+dx-dis_shift)];
    u_shift = [u, u(:, 1)];
    Z_shift = [Z, Z(:, 1)];
else
    %X_shift = [X_shift, X_shift(X_shift<=0) + L + dis_shift];  %enforce periodicity, but need unique points so keep both ends
    X_shift = [X_shift, X_shift(:, inds_outofbox) + L + dis_shift];  %enforce periodicity, but need unique points so keep both ends
    u_shift = [u, u(:, inds_outofbox)];
    Z_shift = [Z, Z(:, inds_outofbox)];
end



%can use gridded interpolant
u_interp = interp2(X_shift, Z_shift, u_shift, X, Z_newgrid);

end



