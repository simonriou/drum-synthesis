% Parameters
Fs = 44100; % sample rate
T = 2; % seconds
T60 = 8; % seconds

NT = round(T * Fs); % time steps
k = 1/Fs; % time step
lambda = 1/sqrt(2); % courant number
gamma = 400;
sigma0 = 6*log(10)/T60;

h = gamma*k/lambda; % grid spacing
epsilon = 1.3; % aspect ratio
Nx = floor(sqrt(epsilon)/h);
Ny = floor(1/(sqrt(epsilon)*h));
h = sqrt(epsilon)/Nx; lambda = gamma*k/h;

Lx = sqrt(epsilon); Ly = 1/sqrt(epsilon);

N = floor(Nx/3);           % Number of modes in x
M = floor(Ny/3);           % Number of modes in y

% Spatial grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[xx, yy] = meshgrid(x, y);

% Initial displacement f(x, y) as raised cosine bump
x0 = Lx/2;
y0 = Ly/2;
R = 0.2;
f = @(x, y) ((x - x0).^2 + (y - y0).^2 <= R^2) .* ...
           0.5 .* (1 + cos(pi * sqrt((x - x0).^2 + (y - y0).^2) / R));

F = f(xx, yy);  % Evaluate on grid

% Compute modal coefficients A_{n,m}
A = zeros(N+1, M+1);
for n = 0:N
    for m = 0:M
        cos_n = cos(n*pi*x'/Lx);    % Nx×1
        cos_m = cos(m*pi*y/Ly);     % 1×Ny
        basis = cos_n * cos_m;      % Nx×Ny
        integrand = F .* basis';    % transpose to match
        A(n+1,m+1) = 4/(Lx*Ly) * trapz(y, trapz(x, integrand, 2));
    end
end

% Time loop setup
t_max = T;
dt = k;
frames = floor(t_max / dt);

% Create figure
figure;
colormap(jet);
h_img = imagesc(x, y, zeros(Ny, Nx));
axis equal tight;
colorbar;
xlabel('x'); ylabel('y');
title('2D Damped Wave Equation Solution');

% Precompute modal frequencies with damping
omega_damped = zeros(N+1, M+1);
for n = 0:N
    for m = 0:M
        k2 = (n*pi/Lx)^2 + (m*pi/Ly)^2;
        omega2 = gamma^2 * k2 - sigma0^2;
        if omega2 > 0
            omega_damped(n+1, m+1) = sqrt(omega2);
        else
            omega_damped(n+1, m+1) = 0;  % overdamped
        end
    end
end

% Choose a specific time index
targetIdx = 300;  % for example, time t = 1.2 seconds
t = (targetIdx - 1) * k;
U_xt = zeros(Ny, Nx);
for n = 0:N
    for m = 0:M
        w_nm = omega_damped(n+1, m+1);
        decay = exp(-sigma0 * t);
        osc = cos(w_nm * t);
        mode_shape = cos(n*pi*xx/Lx) .* cos(m*pi*yy/Ly);
        U_xt = U_xt + A(n+1, m+1) * decay * osc * mode_shape;
    end
end

% Plot and save the frame
figure;
colormap(jet);
imagesc(x, y, U_xt);
axis equal tight;
colorbar;
xlabel('x'); ylabel('y');
title(sprintf('Frame at t = %.3f seconds', t));

% Save to file
filename = sprintf('modal_wave_frame_t%.3fs.png', t);
exportgraphics(gca, filename);  % or use saveas(gcf, filename)
disp(['Saved frame to ', filename]);

% Animation loop
for idx = 1:frames
    t = (idx-1)*dt;
    U_xt = zeros(Ny, Nx);
    for n = 0:N
        for m = 0:M
            w_nm = omega_damped(n+1, m+1);
            decay = exp(-sigma0 * t);
            osc = cos(w_nm * t);
            mode_shape = cos(n*pi*xx/Lx) .* cos(m*pi*yy/Ly);
            U_xt = U_xt + A(n+1, m+1) * decay * osc * mode_shape;
        end
    end
    set(h_img, 'CData', U_xt);
    title(sprintf('2D Damped Wave Equation at time step: %.3f', idx));
    drawnow;
end