clearvars
clc
set(0,'defaulttextInterpreter','latex')
mfile_name          = mfilename('fullpath');
if contains(mfile_name,'LiveEditorEvaluationHelper')
    mfile_name=matlab.desktop.editor.getActiveFilename;
end
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
cd('..')
addpath(genpath(pwd)) 

% --- Physical parameters
R0 = 2439.36;         % Radius in km
rho0 = 5429.3;        % Density in kg/m^3
mu0 = 55e9;       % Shear modulus in Pa
eta0 = 3.16e18;       % Viscosity in Pa.s
G0 = 6.67430e-11;     % Gravitational constant

% --- Define ranges
mu_range = logspace(-2, 0, 20) * mu0;            % Shear modulus range
eta_range = logspace(-5, 0, 20) * eta0;            % Viscosity range

% --- Initialize arrays for k2 Love numbers (real and imaginary)
k2_real = zeros(length(mu_range), length(eta_range));
k2_imag = zeros(length(mu_range), length(eta_range));

% Orbital parameters
T = 87.969216879 * 24 * 3600; % Period in seconds
omega0 = 2 * pi / T;

% Loop over mu and eta
for i = 1:length(mu_range)
    for j = 1:length(eta_range)
        % Define the interior model
        Interior_Model(1).R0 = 1e-3 * R0; % Dummy inner layer
        Interior_Model(1).rho0 = rho0;

        Interior_Model(2).R0 = R0;
        Interior_Model(2).rho0 = rho0;
        Interior_Model(2).mu0 = mu_range(i);
        Interior_Model(2).eta0 = eta_range(j);

        % Forcing
        Forcing(1).Td = T;
        Forcing(1).n = 2;
        Forcing(1).m = 0;
        Forcing(1).F = 1;

        % Numerics setup
        Numerics.Nlayers = length(Interior_Model);
        Numerics.method = 'variable';
        Numerics.Nrbase = 200;
        Numerics.parallel_sol = 0;
        Numerics.parallel_gen = 0;
        Numerics.perturbation_order = 2;
        Numerics.solution_cutoff = 12;
        Numerics.load_couplings = 1;
        Numerics.Nenergy = 12;
        Numerics.rheology_cutoff = 2;

        % Setup boundary and rheology
        [Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model, 'verbose');
        Interior_Model_aux = get_rheology(Interior_Model, Numerics, Forcing);

        % Compute tidal response
        [Love_Spectra,~] = get_Love(Interior_Model_aux, Forcing, Numerics, 'verbose');

        % Store results
        k2_real(i,j) = real(Love_Spectra.k);
        k2_imag(i,j) = imag(Love_Spectra.k);

        % fprintf('mu=%.3e, eta=%.3e, k2=%.4f + %.4fi\n', mu_range(i), eta_range(j), k2_real(i,j), k2_imag(i,j));
    end
end
%%
% --- Plotting ---
figure;
subplot(1,2,1)
pcolor(eta_range / eta0, mu_range / mu0, k2_real)
contourf(eta_range / eta0, mu_range / mu0, k2_real)
shading interp
set(gca, 'XScale', 'log', 'Fontsize', 14)
set(gca, 'YScale', 'log', 'Fontsize', 14)
xlabel('Viscosity $\eta / \eta_0$')
ylabel('Shear modulus $\mu / \mu_0$')
title('Real part of $k_2$', 'FontSize',14)
colorbar

subplot(1,2,2)
pcolor(eta_range / eta0, mu_range / mu0, k2_imag)
contourf(eta_range / eta0, mu_range / mu0, k2_imag)
shading interp
set(gca, 'XScale', 'log', 'Fontsize', 14)
set(gca, 'YScale', 'log', 'Fontsize', 14)
xlabel('Viscosity $\eta / \eta_0$')
ylabel('Shear modulus $\mu / \mu_0$')
title('Imaginary part of $k_2$', 'FontSize', 14)
colorbar
