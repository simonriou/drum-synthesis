%% Parameters
Fs = 44100; % sample rate
T = 2; % seconds
T60 = 8; % seconds

NT = round(T * Fs); % time steps
k = 1/Fs; % time step
lambda = 1/sqrt(2); % courant nb
gamma = 400;
sigma0 = 6*log(10)/T60;

h = gamma*k/lambda; % find grid spacing
epsilon = 1; % aspect ratio
Nx = floor(sqrt(epsilon)/h);
Ny = floor(1/(sqrt(epsilon)*h));
h = sqrt(epsilon)/Nx; lambda = gamma*k/h;

Lx = sqrt(epsilon); Ly = 1/sqrt(epsilon);

%% Listening point & interpolation
x_probe = 0.3;
y_probe = 0.4;
ix = floor(x_probe / h); iy = floor(y_probe / h);
ix = min(max(ix, 1), Nx); iy = min(max(iy, 1), Ny);

x1 = ix * h; y1 = iy * h;
dx = (x_probe - x1) / h; dy = (y_probe - y1) / h;

% Weights
w11 = (1 - dx) * (1 - dy);
w21 = dx * (1 - dy);
w12 = (1 - dx) * dy;
w22 = dx * dy;

% Indices of the 4 corners in flattened vector representation
i11 = iy + (ix-1)*(Ny+1);
i21 = iy + (ix)*(Ny+1);
i12 = (iy+1) + (ix-1)*(Ny+1);
i22 = (iy+1) + (ix)*(Ny+1);

%% Build the Laplacian
function L = build_laplacian(Nx, Ny)

    Ny1 = Ny + 1;
    Nx1 = Nx + 1;
    N = Ny1 * Nx1;
    
    % Preallocate sparse matrix
    L = sparse(N, N);
    
    % Construct diagonal and off-diagonal blocks
    for i = 1:Nx1
        % --- Build diagonal block ---
        d = -4 * ones(Ny1, 1);
        d(1) = -3;
        d(end) = -3;
        if i == 1 || i == Nx1
            d(:) = -3;
            d(1) = -2;
            d(end) = -2;
        end
        main_block = spdiags([ones(Ny1,1), d, ones(Ny1,1)], -1:1, Ny1, Ny1);

        % --- Place diagonal block ---
        row_idx = (i-1)*Ny1 + (1:Ny1);
        col_idx = row_idx;
        L(row_idx, col_idx) = main_block;

        % --- Add identity to upper block (if not last) ---
        if i < Nx1
            next_row_idx = i*Ny1 + (1:Ny1);
            L(row_idx, next_row_idx) = speye(Ny1);
            L(next_row_idx, row_idx) = speye(Ny1);
        end
    end
end

L = 1/h^2 * build_laplacian(Nx, Ny);

%% Initial condition and time integration
[X, Y] = meshgrid(linspace(0, Lx, Nx+1), linspace(0, Ly, Ny+1));

x0 = Lx / 2;
y0 = Ly / 2;
R = 0.2;
R2 = sqrt((X - x0).^2 + (Y - y0).^2);

U0 = 0.5 * (1 + cos(pi * R2 / R));
U0(R2 >= R) = 0;

u_prev = U0(:);
u_curr = u_prev;

audio_out = zeros(1, NT);

n_target = 300;
frame_saved = false;

for n = 1:NT
    u_next = 1/(k*sigma0 + 0.5) * u_curr + (k^2 * gamma^2)/(1 + 2*k*sigma0) * (L * u_curr) + ...
        (2 * sigma0 * k - 1)/(2 * sigma0 * k + 1) * u_prev;

    U = reshape(u_next, Ny+1, Nx+1);
    imagesc(U); colorbar; axis equal tight;
    title(['Time step: ', num2str(n)]);
    drawnow;

    % Save the desired frame
    if n == n_target && ~frame_saved
        filename = sprintf('fdtd_wave_frame_t%.3fs.png', (n-1)*k);
        exportgraphics(gca, filename);  % Save as PNG
        disp(['Saved frame to ', filename]);
        frame_saved = true;
    end

    audio_out(n) = w11*u_curr(i11) + w21*u_curr(i21) + ...
        w12*u_curr(i12) + w22*u_curr(i22);

    u_prev = u_curr;
    u_curr = u_next;
end

soundsc(audio_out, Fs);
plot((0:NT-1)*k, audio_out, 'k'); xlabel('t'); ylabel('u');
title('2D wave equation'); axis tight;