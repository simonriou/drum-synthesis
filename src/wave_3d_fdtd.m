clear; clc;

%% Parameters
Fs = 44100; % sample rate
T = 5; % seconds
T60 = 8; % seconds

NT = round(T * Fs); % time steps
k = 1/Fs; % time step
lambda = 1/sqrt(3); % courant nb
gamma = 340;
sigma0 = 6*log(10)/T60;

h = gamma*k/lambda; % find grid spacing
epsilon = 1; % aspect ratio y to z
Nz = floor(1/h); % we set Lz = 1
Nx = floor(epsilon/h);
Ny = floor(epsilon/h);
h = sqrt(epsilon)/Nx; lambda = gamma*k/h;

Lx = sqrt(epsilon); Ly = 1/sqrt(epsilon); Lz = 1;

%% Listening point & interpolation
x_probe = 0.3;
y_probe = 0.4;
z_probe = 0.5;
ix = floor(x_probe / h); iy = floor(y_probe / h); iz = floor(z_probe / h);
ix = min(max(ix, 1), Nx); iy = min(max(iy, 1), Ny); iz = min(max(iz, 1), Nz);

x1 = ix * h; y1 = iy * h; z1 = iz * h;
dx = (x_probe - x1) / h; dy = (y_probe - y1) / h; dz = (z_probe - z1) / h;

% Weights
w000 = (1 - dx) * (1 - dy) * (1 - dz);
w100 = dx * (1 - dy) * (1 - dz);
w010 = (1 - dx) * dy * (1 - dz);
w110 = dx * dy * (1 - dz);
w001 = (1 - dx) * (1 - dy) * dz;
w101 = dx * (1 - dy) * dz;
w011 = (1 - dx) * dy * dz;
w111 = dx * dy * dz;

% Indices of the 8 corners in flattened vector representation
i000 = iz     + (iy-1)   *(Nz+1) + (ix-1)   *(Ny+1)*(Nz+1);
i100 = iz     + (iy-1)   *(Nz+1) + (ix)     *(Ny+1)*(Nz+1);
i010 = iz     + (iy)     *(Nz+1) + (ix-1)   *(Ny+1)*(Nz+1);
i110 = iz     + (iy)     *(Nz+1) + (ix)     *(Ny+1)*(Nz+1);
i001 = (iz+1) + (iy-1)   *(Nz+1) + (ix-1)   *(Ny+1)*(Nz+1);
i101 = (iz+1) + (iy-1)   *(Nz+1) + (ix)     *(Ny+1)*(Nz+1);
i011 = (iz+1) + (iy)     *(Nz+1) + (ix-1)   *(Ny+1)*(Nz+1);
i111 = (iz+1) + (iy)     *(Nz+1) + (ix)     *(Ny+1)*(Nz+1);

%% Build the Laplacian
function L = build_laplacian(Nx, Ny, Nz)
    Ny1 = Ny + 1;
    Nx1 = Nx + 1;
    Nz1 = Nz + 1;

    Nxy = Nx1 * Ny1;
    N = Nxy * Nz1;
    
    % Build the 2D Laplacian matrix (d2)
    d2 = sparse(Nxy, Nxy);
    for i = 1:Nx1
        % Build diagonal block
        d = -4 * ones(Ny1, 1);
        d(1) = -3;
        d(end) = -3;
        if i == 1 || i == Nx1
            d(:) = -3;
            d(1) = -2;
            d(end) = -2;
        end
        main_block = spdiags([ones(Ny1,1), d, ones(Ny1,1)], -1:1, Ny1, Ny1);

        % Place diagonal block
        row_idx = (i-1)*Ny1 + (1:Ny1);
        col_idx = row_idx;
        d2(row_idx, col_idx) = main_block;

        % Add identity to upper and lower blocks if not the last x layer
        if i < Nx1
            next_row_idx = i*Ny1 + (1:Ny1);
            d2(row_idx, next_row_idx) = speye(Ny1);
            d2(next_row_idx, row_idx) = speye(Ny1);
        end
    end

    % Construct the 3D Laplacian matrix
    x_y_part = kron(speye(Nz1), d2);

    % Diagonal contributions from the z-direction
    z_contrib = zeros(N, 1);
    for k = 1:Nz1
        start = (k-1)*Nxy + 1;
        finish = k*Nxy;
        if k == 1 || k == Nz1
            z_contrib(start:finish) = -1;
        else
            z_contrib(start:finish) = -2;
        end
    end
    Z_diag = spdiags(z_contrib, 0, N, N);

    % Off-diagonal contributions from the z-direction (identity matrices)
    num_upper = Nxy * (Nz1 - 1);
    rows_upper = (1 : num_upper)';
    cols_upper = rows_upper + Nxy;
    vals_upper = ones(num_upper, 1);

    rows_lower = cols_upper;
    cols_lower = rows_upper;
    vals_lower = ones(num_upper, 1);

    Z_offdiag = sparse([rows_upper; rows_lower], [cols_upper; cols_lower], [vals_upper; vals_lower], N, N);

    % Combine all parts to form the 3D Laplacian
    L = x_y_part + Z_diag + Z_offdiag;
end

L = build_laplacian(Nx, Ny, Nz);

%% Initial condition and time integration
[X, Y, Z] = meshgrid(linspace(0, Lx, Nx+1), ...
                    linspace(0, Ly, Ny+1), ...
                    linspace(0, Lz, Nz+1));

x0 = Lx / 2;
y0 = Ly / 2;
z0 = Lz / 2;

R = 0.2;
R3 = sqrt((X - x0).^2 + (Y - y0).^2 + (Z - z0).^2);

U0 = 0.5 * (1 + cos(pi * R3 / R));
U0(R3 >= R) = 0;

u_prev = U0(:);
u_curr = u_prev;

audio_out = zeros(1, NT);

%% Time Loop

% Define z-slice index for visualization
z_slice = round((Nz+1)/2); % center of the volume

% Prepare figure and initialize image object on first iteration
figure;
h_img = [];
h_title = [];

for n = 1:NT
    % Display progress in console
    disp("Iteration " + n + "/" + NT)

    % Time update step
    % u_next = (1 / (1 + sigma0 * k)) * ...
    %      (2 * u_curr - (1 - sigma0 * k) * u_prev + k^2 * gamma^2 * (L * u_curr));
    u_next = 2*u_curr + lambda^2 * L * u_curr - u_prev;

    % Audio output interpolation
    audio_out(n) = ...
        w000*u_curr(i000) + w100*u_curr(i100) + ...
        w010*u_curr(i010) + w110*u_curr(i110) + ...
        w001*u_curr(i001) + w101*u_curr(i101) + ...
        w011*u_curr(i011) + w111*u_curr(i111);

    % Efficient plotting every 100 time steps
    % if mod(n, 100) == 0
    %     U = reshape(u_next, Nz+1, Ny+1, Nx+1); % Proper reshaping
    %     U_slice = squeeze(U(z_slice, :, :));  % Extract middle z-slice

    %     if isempty(h_img)
    %         h_img = imagesc(U_slice);
    %         colorbar;
    %         axis equal tight;
    %         h_title = title(['Time step: ', num2str(n), ' | z-slice: ', num2str(z_slice)]);
    %     else
    %         set(h_img, 'CData', U_slice);
    %         set(h_title, 'String', ['Time step: ', num2str(n), ' | z-slice: ', num2str(z_slice)]);
    %     end

    %     drawnow limitrate
    % end

    % Advance time step
    u_prev = u_curr;
    u_curr = u_next;
end

soundsc(audio_out, Fs);
plot((0:NT-1)*k, audio_out, 'k'); xlabel('t'); ylabel('u');
title('3D wave equation'); axis tight;