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

% Define constant radii
R_outer_const = 2010;
R_mantle_const = 2410;
R_crust_const = 2439.36;
R_inner_const = 1180;

%Define constant density
rho_inner_const = 7380;
rho_outer_const = 6470;
rho_mantle_const = 3100;
rho_crust_const = 2850;

% Define constant shear moduli 
mu_mantle_const = 6e10;
mu_crust_const = 5.5e10;
mu_inner_const = 10e10;

% Define constant viscosity
eta_mantle_const = 1e18;
eta_crust_const = 1e20;
eta_inner_const = 1e23;


% Define radii range
N = 21;
x_0 = 0.9;
x_n = 1.1;
R_inner_range = linspace(0.9, 1.1, N) * 1180;
R_outer_range = linspace(0.9, 1.1, N) * 2010;
R_mantle_range = linspace(0.9, 1.1, N) * 2410;

% Define density change
rho_inner_range = linspace(0.9, 1.1, N)*rho_inner_const;
rho_outer_range = linspace(0.9, 1.1, N)*rho_outer_const;
rho_mantle_range = linspace(0.9, 1.1, N)*rho_mantle_const;
rho_crust_range = linspace(0.9, 1.1, N)*rho_crust_const;

% Define shear moduli range
mu_inner_range = linspace(0.9, 1.1, N) * 10e10;
mu_mantle_range = linspace(0.9, 1.1, N) * 6e10;
mu_crust_range = linspace(0.9, 1.1, N) * 5.5e10;

% Define viscocity range
% eta_inner_range = logspace(-5, 5, N) * eta_inner_const;
% eta_mantle_range = logspace(-5, 5, N) * eta_mantle_const;
% eta_crust_range = logspace(-5, 5, N) * eta_crust_const;

% Define viscocity range
eta_inner_range = linspace(0.9, 1.1, N) * eta_inner_const;
eta_mantle_range = linspace(0.9, 1.1, N) * eta_mantle_const;
eta_crust_range = linspace(0.9, 1.1, N) * eta_crust_const;

% Initialize result arrays
k2_R_inner = zeros(N,1);
k2_R_outer = zeros(N,1);
k2_R_mantle = zeros(N,1);
k2_rho_inner = zeros(N,1);
k2_rho_outer = zeros(N,1);
k2_rho_mantle = zeros(N,1);
k2_rho_crust = zeros(N,1);
k2_mu_inner = zeros(N,1);
k2_mu_crust = zeros(N,1);
k2_mu_mantle = zeros(N,1);
k2_eta_inner = zeros(N,1);
k2_eta_crust = zeros(N,1);
k2_eta_mantle = zeros(N,1);
%% LOOP 1: Vary inner core radius
for i = 1:N
    Interior_Model = setup_radii(R_inner_range(i), R_outer_const, R_mantle_const, R_crust_const);
    [k2_R_inner(i), h2_R_inner(i)] = compute_k2(Interior_Model);
end

%% LOOP 2: Vary outer core radius
for i = 1:N
    Interior_Model = setup_radii(R_inner_const, R_outer_range(i), R_mantle_const, R_crust_const);
    [k2_R_outer(i), h2_R_outer(i)] = compute_k2(Interior_Model);
end

%% LOOP 3: Vary mantle radius
for i = 1:N
    Interior_Model = setup_radii(R_inner_const, R_outer_const, R_mantle_range(i), R_crust_const);
    [k2_R_mantle(i), h2_R_mantle(i)] = compute_k2(Interior_Model);
end

%% LOOP 4: Vary inner core shear modulus
for i = 1:N
    Interior_Model = setup_shear_moduli(mu_inner_range(i), mu_mantle_const, mu_crust_const);
    [k2_mu_inner(i), h2_mu_inner(i)] = compute_k2(Interior_Model);
end

%% LOOP 5: Vary outer core mantle shear modulus
for i = 1:N
    Interior_Model = setup_shear_moduli(mu_inner_const, mu_mantle_range(i), mu_crust_const);
    [k2_mu_mantle(i), h2_mu_mantle(i)] = compute_k2(Interior_Model);
end

%% LOOP 6: Vary crust shear modulus
for i = 1:N
    Interior_Model = setup_shear_moduli(mu_inner_const, mu_mantle_const, mu_crust_range(i));
    [k2_mu_crust(i), h2_mu_crust(i)] = compute_k2(Interior_Model);
end

%% LOOP 7: Vary inner core viscosity
for i = 1:N
    Interior_Model = setup_viscosities(eta_inner_range(i), eta_mantle_const, eta_crust_const);
    [k2_eta_inner(i), h2_eta_inner(i)] = compute_k2(Interior_Model);
end

%% LOOP 8: Vary mantle viscosity
for i = 1:N
    Interior_Model = setup_viscosities(eta_inner_const, eta_mantle_range(i), eta_crust_const);
    [k2_eta_mantle(i), h2_eta_mantle(i)] = compute_k2(Interior_Model);
end

%% LOOP 9: Vary crust viscosity
for i = 1:N
    Interior_Model = setup_viscosities(eta_inner_const, eta_mantle_const, eta_crust_range(i));
    [k2_eta_crust(i), h2_eta_crust(i)] = compute_k2(Interior_Model);
end

%% LOOP 10: Vary inner core density
for i = 1:N
    Interior_Model = setup_densities(rho_inner_range(i), rho_outer_const, rho_mantle_const, rho_crust_const);
    [k2_rho_inner(i), h2_rho_inner(i)] = compute_k2(Interior_Model);
end

%% LOOP 11: Vary outer core density
for i = 1:N
    Interior_Model = setup_densities(rho_inner_const, rho_outer_range(i), rho_mantle_const, rho_crust_const);
    [k2_rho_outer(i), h2_rho_outer(i)] = compute_k2(Interior_Model);
end

%% LOOP 12: Vary mantle density
for i = 1:N
    Interior_Model = setup_densities(rho_inner_const, rho_outer_const, rho_mantle_range(i), rho_crust_const);
    [k2_rho_mantle(i), h2_rho_mantle(i)] = compute_k2(Interior_Model);
end

%% LOOP 13: Vary crust density
for i = 1:N
    Interior_Model = setup_densities(rho_inner_const, rho_outer_const, rho_mantle_const, rho_crust_range(i));
    [k2_rho_crust(i), h2_rho_crust(i)] = compute_k2(Interior_Model);
end

%% Calculating derivatives
h = (x_n-x_0)/N;
dk2dR_inner = gradient(k2_R_inner, h);
dk2dR_outer = gradient(k2_R_outer, h);
dk2dR_mantle = gradient(k2_R_mantle, h);
dk2dmu_inner = gradient(k2_mu_inner, h);
dk2dmu_mantle = gradient(k2_mu_mantle, h);
dk2dmu_crust = gradient(k2_mu_crust, h);
dk2deta_inner = gradient(k2_eta_inner, h);
dk2deta_mantle = gradient(k2_eta_mantle, h);
dk2deta_crust = gradient(k2_eta_crust, h);
dk2drho_inner = gradient(k2_rho_inner, h);
dk2drho_outer = gradient(k2_rho_outer, h);
dk2drho_mantle = gradient(k2_rho_mantle, h);
dk2drho_crust = gradient(k2_rho_crust, h);

fprintf('dk2dR_inner:%.3e dk2dR_outer:%.3e dk2dR_mantle:%.3e dk2drho_inner:%.3e d dk2drho_outer:%.3e d dk2drho_mantle:%.3e d dk2drho_crust:%.3e d k2dmu_inner:%.3e dk2dmu_mantle:%.3e dk2dmu_crust:%.3e dk2deta_inner:%.3e dk2deta_mantle:%.3e dk2deta_crust:%.3e\n', dk2dR_inner(11), dk2dR_outer(11), dk2dR_mantle(11), dk2drho_inner(11), dk2drho_outer(11), dk2drho_mantle(11), dk2drho_crust(11), dk2dmu_inner(11), dk2dmu_mantle(11), dk2dmu_crust(11), dk2deta_inner(11), dk2deta_mantle(11), dk2deta_crust(11))

%% Plotting k2
close all;
figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(R_inner_range / R_inner_const, real(k2_R_inner), 'r-', 'LineWidth', 2)
plot(R_outer_range / R_outer_const, real(k2_R_outer), 'r--', 'LineWidth', 2)
plot(R_mantle_range / R_mantle_const, real(k2_R_mantle), 'r.', 'LineWidth', 2)
plot(rho_inner_range / rho_inner_const, real(k2_rho_inner), 'c-', 'Linewidth', 2)
plot(rho_outer_range / rho_outer_const, real(k2_rho_outer), 'c--', 'Linewidth', 2)
plot(rho_mantle_range / rho_mantle_const, real(k2_rho_mantle), 'c.', 'Linewidth', 2)
plot(rho_crust_range / rho_crust_const, real(k2_rho_crust), 'c-.', 'Linewidth', 2)
plot(mu_inner_range / mu_inner_const, real(k2_mu_inner), 'g-', 'LineWidth', 2)
plot(mu_mantle_range / mu_mantle_const, real(k2_mu_mantle), 'g--', 'LineWidth', 2)
plot(mu_crust_range / mu_crust_const, real(k2_mu_crust), 'g.', 'LineWidth', 2)
plot(eta_inner_range / eta_inner_const, real(k2_eta_inner), 'b-', 'LineWidth', 2)
plot(eta_mantle_range / eta_mantle_const, real(k2_eta_mantle), 'b--', 'LineWidth', 2)
plot(eta_crust_range / eta_crust_const, real(k2_eta_crust), 'b.', 'LineWidth', 2)
xlabel('Normalized Parameter $x / x_0$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0.9,1.1]);
ylim([0, 1]);
ylabel('Real part of $k_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Radius Inner Core', 'Radius Outer Core', 'Radius Mantle', 'Density Inner Core', 'Density Outer Core', 'Density Mantle', 'Density Crust','Shear Modulus Inner Core', 'Shear Modulus Mantle', 'Shear Modulus Crust', 'Viscosity Inner Core', 'Viscosity Mantle', 'Viscosity Crust'}, 'Location', 'best')
grid on
savefig('Assignment2/sensitivity.fig')
saveas(gcf, 'Assignment2/sensitivity.svg')

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(R_inner_range / R_inner_const, imag(k2_R_inner), 'r-', 'LineWidth', 2)
plot(R_outer_range / R_outer_const, imag(k2_R_outer), 'r--', 'LineWidth', 2)
plot(R_mantle_range / R_mantle_const, imag(k2_R_mantle), 'r.', 'LineWidth', 2)
plot(rho_inner_range / rho_inner_const, imag(k2_rho_inner), 'c-', 'Linewidth', 2)
plot(rho_outer_range / rho_outer_const, imag(k2_rho_outer), 'c--', 'Linewidth', 2)
plot(rho_mantle_range / rho_mantle_const, imag(k2_rho_mantle), 'c.', 'Linewidth', 2)
plot(rho_crust_range / rho_crust_const, imag(k2_rho_crust), 'c-.', 'Linewidth', 2)
plot(mu_inner_range / mu_inner_const, imag(k2_mu_inner), 'g-', 'LineWidth', 2)
plot(mu_mantle_range / mu_mantle_const, imag(k2_mu_mantle), 'g--', 'LineWidth', 2)
plot(mu_crust_range / mu_crust_const, imag(k2_mu_crust), 'g.', 'LineWidth', 2)
plot(eta_inner_range / eta_inner_const, imag(k2_eta_inner), 'b-', 'LineWidth', 2)
plot(eta_mantle_range / eta_mantle_const, imag(k2_eta_mantle), 'b--', 'LineWidth', 2)
plot(eta_crust_range / eta_crust_const, imag(k2_eta_crust), 'b.', 'LineWidth', 2)
xlabel('Normalized Parameter $x / x_0$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0.9,1.1]);
% ylim([0, 1]);
ylabel('Imaginary part of $k_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Radius Inner Core', 'Radius Outer Core', 'Radius Mantle', 'Density Inner Core', 'Density Outer Core', 'Density Mantle', 'Density Crust','Shear Modulus Inner Core', 'Shear Modulus Mantle', 'Shear Modulus Crust', 'Viscosity Inner Core', 'Viscosity Mantle', 'Viscosity Crust'}, 'Location', 'best')
grid on
savefig('Assignment2/sensitivity_imag.fig')
saveas(gcf, 'Assignment2/sensitivity_imag.svg')

%% Plotting h2
figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(R_inner_range / R_inner_const, real(h2_R_inner), 'r-', 'LineWidth', 2)
plot(R_outer_range / R_outer_const, real(h2_R_outer), 'r--', 'LineWidth', 2)
plot(R_mantle_range / R_mantle_const, real(h2_R_mantle), 'r.', 'LineWidth', 2)
plot(rho_inner_range / rho_inner_const, real(h2_rho_inner), 'c-', 'Linewidth', 2)
plot(rho_outer_range / rho_outer_const, real(h2_rho_outer), 'c--', 'Linewidth', 2)
plot(rho_mantle_range / rho_mantle_const, real(h2_rho_mantle), 'c.', 'Linewidth', 2)
plot(rho_crust_range / rho_crust_const, real(h2_rho_crust), 'c-.', 'Linewidth', 2)
plot(mu_inner_range / mu_inner_const, real(h2_mu_inner), 'g-', 'LineWidth', 2)
plot(mu_mantle_range / mu_mantle_const, real(h2_mu_mantle), 'g--', 'LineWidth', 2)
plot(mu_crust_range / mu_crust_const, real(h2_mu_crust), 'g.', 'LineWidth', 2)
plot(eta_inner_range / eta_inner_const, real(h2_eta_inner), 'b-', 'LineWidth', 2)
plot(eta_mantle_range / eta_mantle_const, real(h2_eta_mantle), 'b--', 'LineWidth', 2)
plot(eta_crust_range / eta_crust_const, real(h2_eta_crust), 'b.', 'LineWidth', 2)
xlabel('Normalized Parameter $x / x_0$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0.9,1.1]);
% ylim([0, 1]);
ylabel('Real part of $h_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Radius Inner Core', 'Radius Outer Core', 'Radius Mantle', 'Density Inner Core', 'Density Outer Core', 'Density Mantle', 'Density Crust','Shear Modulus Inner Core', 'Shear Modulus Mantle', 'Shear Modulus Crust', 'Viscosity Inner Core', 'Viscosity Mantle', 'Viscosity Crust'}, 'Location', 'best')
grid on
savefig('Assignment2/sensitivity_h.fig')
saveas(gcf, 'Assignment2/sensitivity_h.svg')

figure('units', 'points', 'Position', [0, 0, 455.2441, 0.8*455.2441]);
hold on
% Normalized x-axis: divide each by its central value
plot(R_inner_range / R_inner_const, imag(h2_R_inner), 'r-', 'LineWidth', 2)
plot(R_outer_range / R_outer_const, imag(h2_R_outer), 'r--', 'LineWidth', 2)
plot(R_mantle_range / R_mantle_const, imag(h2_R_mantle), 'r.', 'LineWidth', 2)
plot(rho_inner_range / rho_inner_const, imag(h2_rho_inner), 'c-', 'Linewidth', 2)
plot(rho_outer_range / rho_outer_const, imag(h2_rho_outer), 'c--', 'Linewidth', 2)
plot(rho_mantle_range / rho_mantle_const, imag(h2_rho_mantle), 'c.', 'Linewidth', 2)
plot(rho_crust_range / rho_crust_const, imag(h2_rho_crust), 'c-.', 'Linewidth', 2)
plot(mu_inner_range / mu_inner_const, imag(h2_mu_inner), 'g-', 'LineWidth', 2)
plot(mu_mantle_range / mu_mantle_const, imag(h2_mu_mantle), 'g--', 'LineWidth', 2)
plot(mu_crust_range / mu_crust_const, imag(h2_mu_crust), 'g.', 'LineWidth', 2)
plot(eta_inner_range / eta_inner_const, imag(h2_eta_inner), 'b-', 'LineWidth', 2)
plot(eta_mantle_range / eta_mantle_const, imag(h2_eta_mantle), 'b--', 'LineWidth', 2)
plot(eta_crust_range / eta_crust_const, imag(h2_eta_crust), 'b.', 'LineWidth', 2)
xlabel('Normalized Parameter $x / x_0$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0.9,1.1]);
% ylim([0, 1]);
ylabel('Imaginary part of $h_2$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'Radius Inner Core', 'Radius Outer Core', 'Radius Mantle', 'Density Inner Core', 'Density Outer Core', 'Density Mantle', 'Density Crust','Shear Modulus Inner Core', 'Shear Modulus Mantle', 'Shear Modulus Crust', 'Viscosity Inner Core', 'Viscosity Mantle', 'Viscosity Crust'}, 'Location', 'best')
grid on
savefig('Assignment2/sensitivity_h_imag.fig')
saveas(gcf, 'Assignment2/sensitivity_h_imag.svg')

%% Functions
function Interior_Model = setup_radii(R_inner, R_outer, R_mantle, R_crust)
    % Inner core (solid)
    Interior_Model(1).R0 = 1e-3 * R_inner; % dummy layer
    Interior_Model(1).rho0 = 7380;

    Interior_Model(2).R0 = R_inner;
    Interior_Model(2).rho0 = 7380;
    Interior_Model(2).mu0 = 10.0e10;
    Interior_Model(2).eta0 = 1e20;

    % Outer core
    Interior_Model(3).R0 = R_outer;
    Interior_Model(3).rho0 = 6470;
    Interior_Model(3).ocean = 1;

    % Mantle
    Interior_Model(4).R0 = R_mantle;
    Interior_Model(4).rho0 = 3100;
    Interior_Model(4).mu0 = 6e10;
    Interior_Model(4).eta0 = 1e18;

    % Crust
    Interior_Model(5).R0 = R_crust;
    Interior_Model(5).rho0 = 2850;
    Interior_Model(5).mu0 = 5.5e10;
    Interior_Model(5).eta0 = 1e23;
end

function Interior_Model = setup_densities(rho_inner, rho_outer, rho_mantle, rho_crust)
    % Inner core (solid)
    Interior_Model(1).R0 = 1e-3 * 1180; % dummy layer
    Interior_Model(1).rho0 = 7380;

    Interior_Model(2).R0 = 1180;
    Interior_Model(2).rho0 = rho_inner;
    Interior_Model(2).mu0 = 10.0e10;
    Interior_Model(2).eta0 = 1e20;

    % Outer core
    Interior_Model(3).R0 = 2010;
    Interior_Model(3).rho0 = rho_outer;
    Interior_Model(3).ocean = 1;

    % Mantle
    Interior_Model(4).R0 = 2410;
    Interior_Model(4).rho0 = rho_mantle;
    Interior_Model(4).mu0 = 6e10;
    Interior_Model(4).eta0 = 1e18;

    % Crust
    Interior_Model(5).R0 = 2439.36;
    Interior_Model(5).rho0 = rho_crust;
    Interior_Model(5).mu0 = 5.5e10;
    Interior_Model(5).eta0 = 1e23;
end

function Interior_Model = setup_shear_moduli(mu_inner, mu_mantle, mu_crust)
    % Inner core (solid)
    Interior_Model(1).R0 = 1e-3 * 1180; % dummy layer
    Interior_Model(1).rho0 = 7380;

    Interior_Model(2).R0 = 1180;
    Interior_Model(2).rho0 = 7380;
    Interior_Model(2).mu0 = mu_inner;
    Interior_Model(2).eta0 = 1e20;

    % Outer core
    Interior_Model(3).R0 = 2010;
    Interior_Model(3).rho0 = 6470;
    Interior_Model(3).ocean = 1;

    % Mantle
    Interior_Model(4).R0 = 2410;
    Interior_Model(4).rho0 = 3100;
    Interior_Model(4).mu0 = mu_mantle;
    Interior_Model(4).eta0 = 1e18;

    % Crust
    Interior_Model(5).R0 = 2439.36;
    Interior_Model(5).rho0 = 2850;
    Interior_Model(5).mu0 = mu_crust;
    Interior_Model(5).eta0 = 1e23;
end

function Interior_Model = setup_viscosities(eta_inner, eta_mantle, eta_crust)
    % Inner core (solid)
    Interior_Model(1).R0 = 1e-3 * 1180; % dummy layer
    Interior_Model(1).rho0 = 7380;

    Interior_Model(2).R0 = 1180;
    Interior_Model(2).rho0 = 7380;
    Interior_Model(2).mu0 = 10e10;
    Interior_Model(2).eta0 = eta_inner;

    % Outer core
    Interior_Model(3).R0 = 2010;
    Interior_Model(3).rho0 = 6470;
    Interior_Model(3).ocean = 1;

    % Mantle
    Interior_Model(4).R0 = 2410;
    Interior_Model(4).rho0 = 3100;
    Interior_Model(4).mu0 = 6e10;
    Interior_Model(4).eta0 = eta_mantle;

    % Crust
    Interior_Model(5).R0 = 2439.36;
    Interior_Model(5).rho0 = 2850;
    Interior_Model(5).mu0 = 5.5e10;
    Interior_Model(5).eta0 = eta_crust;
end

function [k2_value, h2_value] = compute_k2(Interior_Model) % does in fact campoute k2 and h2
    T = 87.969216879*24*3600;
    Forcing.Td = T;
    Forcing.n = 2; 
    Forcing.m = 0; 
    Forcing.F = 1;

    Numerics.Nlayers = length(Interior_Model);
    Numerics.method = 'variable';
    Numerics.Nrbase = 200;
    Numerics.parallel_sol = 0;
    Numerics.parallel_gen = 0;
    Numerics.perturbation_order = 2;

    [Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model, 'verbose');
    Interior_Model_aux = get_rheology(Interior_Model, Numerics, Forcing);
    [Love_Spectra, ~] = get_Love(Interior_Model_aux, Forcing, Numerics, 'verbose');
    k2_value = Love_Spectra.k;
    h2_value = Love_Spectra.h;
end