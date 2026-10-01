function u_shifted = timeshift_soln(c, dt, L, numshift, X, Z, u, Z_newgrid)
%shift quantity u back c*dt numshift times on periodic grid of length L

%%%% CLAUDE
    % TIMESHIFT_SOLN - Shift solution using spectral (FFT) method
    % Shift quantity u by distance c*dt*numshift on periodic grid of length L
    %
    % This version uses FFT for exact periodic shifting (no interpolation error)
    
    % Total shift distance
    dis_shift = numshift * c * dt;

    % Get dimensions
    [r, col] = size(u);

    % Handle 1D case (row or column vector)
    if r == 1 || col == 1
        % 1D shift
        N = length(u);
        dx = X(1, 2) - X(1, 1);

        % Ensure u is a row vector for FFT
        is_column = (size(u, 1) > 1);
        if is_column
            u = u.';
        end

        % Wavenumbers
        if mod(N, 2) == 0
            k = [0:N/2-1, 0, -N/2+1:-1] * (2*pi/L);
        else
            k = [0:(N-1)/2, -(N-1)/2:-1] * (2*pi/L);
        end

        % FFT shift
        u_hat = fft(u);
        phase_shift = exp(-1i * k * dis_shift);
        u_shifted = real(ifft(u_hat .* phase_shift));

        % Return in original orientation
        if is_column
            u_shifted = u_shifted.';
        end

    else
        % 2D shift (only in x-direction, assuming X varies along columns)
        % Shift each row independently

        dx = X(1, 2) - X(1, 1);
        N = col;  % Number of points in x-direction

        % Wavenumbers
        if mod(N, 2) == 0
            k = [0:N/2-1, 0, -N/2+1:-1] * (2*pi/L);
        else
            k = [0:(N-1)/2, -(N-1)/2:-1] * (2*pi/L);
        end

        % Initialize output
        u_shifted = zeros(size(u));

        % Shift each row
        for i = 1:r
            u_row = u(i, :);
            u_hat = fft(u_row);
            phase_shift = exp(-1i * k * dis_shift);
            u_shifted(i, :) = real(ifft(u_hat .* phase_shift));
        end
    end
end


%%%% MINE 
% dis_shift = numshift*c*dt;
% [r, col] = size(X);
% dx = X(1, 2)-X(1, 1);
% 
% X_shift = X - dis_shift; %shift back how far wave propogated in numshift time steps
% [~, inds_outofbox] = find(X_shift<=0);
% [~, inds_outofbox_right] = find(X_shift>=max(max(X)));
% inds_outofbox=unique(inds_outofbox);
% inds_outofbox_right=unique(inds_outofbox_right);
% 
% if isempty(inds_outofbox) & isempty(inds_outofbox_right)
%     X_shift = [X_shift, ones(r, 1)*(L+dx-dis_shift)];
%     u_shift = [u, u(:, 1)];
%     Z_shift = [Z, Z(:, 1)];
% elseif isempty(inds_outofbox_right)
%     %X_shift = [X_shift, X_shift(X_shift<=0) + L + dis_shift];  %enforce periodicity, but need unique points so keep both ends
%     X_shift = [X_shift, X_shift(:, inds_outofbox) + L + dis_shift];  %enforce periodicity, but need unique points so keep both ends
%     u_shift = [u, u(:, inds_outofbox)];
%     Z_shift = [Z, Z(:, inds_outofbox)];
% else
%      X_shift = [X_shift(:, inds_outofbox_right) - L + dis_shift, X_shift];  %enforce periodicity, but need unique points so keep both ends
%     u_shift = [u(:, inds_outofbox_right), u];
%     if(r>1 && col>1)
%         Z_shift = [Z(:, inds_outofbox_right), Z];  
%     end
% end
% 
% 
% if(r>1 && col>1)
%     %can use gridded interpolant
%     u_interp = interp2(X_shift, Z_shift, u_shift, X, Z_newgrid);
% else
%     u_interp = interp1(X_shift, u_shift, X, "spline");
% end
% 
% u_shifted = u_interp;

