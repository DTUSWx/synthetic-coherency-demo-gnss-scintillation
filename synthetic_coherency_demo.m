%% synthetic_coherency_demo.m
% Synthetic demonstration for the coherency-based GNSS scintillation paper.
%
% Purpose
% -------
% This script generates a synthetic dual-polarization dataset that
% demonstrates three key claims of the manuscript:
%
% 1) The same underlying propagated coherency state can produce different
%    scalar receiver observables when measured through different receiver
%    response matrices.
% 2) Receiver de-embedding can recover receiver-de-embedded coherency
%    quantities that agree across receivers.
% 3) Basis-invariant quantities such as degree of polarization (DoP) and
%    polarization entropy remain constant under unitary polarization
%    rotation, but change under effective mixing / depolarization.
%
% Regimes
% -------
% Regime I   : baseline unitary trajectory
% Regime II  : still unitary, but through a more receiver-sensitive
%              polarization region
% Regime III : mixing / depolarization
%
% Output
% ------
% A four-panel figure suitable as Figure 4, plus publication-quality
% exports (PDF, PNG, EPS, SVG if supported, TIFF if supported, and MATLAB FIG)
% in a subfolder located in the same folder as this script.
%
% Notes
% -----
% - The dataset is synthetic but physically motivated.
% - The field is modeled at the coherency level rather than by full GNSS
%   tracking-loop simulation.
% - The script is self-contained and requires no toolboxes.
%
% Publication note
% ----------------
% This figure is intended as a compact conceptual demonstration of
% receiver dependence, de-embedding recovery, and invariant behavior.
% It does not replace observational validation with calibrated
% dual-polarization GNSS scintillation data.
%
% -------------------------------------------------------------------------

clear; close all; clc;

%% -----------------------------
% User-adjustable settings
% -----------------------------
N = 1200;                  % number of epochs
dt = 0.1;                  % seconds per epoch
t = (0:N-1)' * dt;         % time vector
winLen = 120;              % window length in samples for local scalar metrics
smoothWin = 9;             % odd-length moving average for display smoothing

% Regime boundaries
idx1 = 1:400;              % Regime I
idx2 = 401:800;            % Regime II
idx3 = 801:1200;           % Regime III

% Base total intensity
S0_base = 1.0;

% Initial degree of polarization in the unitary regimes
DoP0 = 0.85;

% Receiver noise power
noisePowA = 0.0035;
noisePowB = 0.0035;

rng(42); % reproducibility

%% -----------------------------
% Plot styling defaults
% -----------------------------
% These values were increased so that Figure 4 better matches the apparent
% font scale used in Figure 5 and remains legible after journal reduction.
fontName = 'Helvetica';

lwMain  = 2.6;
lwRef   = 2.3;
lwAxis  = 1.2;
lwXline = 1.2;

fsTick = 17;   % tick labels
fsAx   = 19;   % axis labels
fsTit  = 20;   % panel titles
fsLeg  = 15;   % legends
fsAnn  = 15;   % in-panel annotations

mkSmall = 8.0;
mkMed   = 8.8;

%% -----------------------------
% Pauli matrices and helpers
% -----------------------------
I2 = eye(2);
sigma1 = [0 1; 1 0];
sigma2 = [0 -1i; 1i 0];
sigma3 = [1 0; 0 -1];

% Hermitian symmetrization helper
herm = @(X) 0.5 * (X + X');

% Simple moving average helper for plotting
movavg = @(x,w) conv(x, ones(w,1)/w, 'same');

% Build coherency from Stokes parameters
J_from_stokes = @(S0,S) 0.5 * (S0*I2 + S(1)*sigma1 + S(2)*sigma2 + S(3)*sigma3);

% Compute Stokes from coherency
stokes_from_J = @(J) [ ...
    real(trace(J)), ...
    real(trace(J*sigma1)), ...
    real(trace(J*sigma2)), ...
    real(trace(J*sigma3)) ];

% Safe DoP from coherency
dop_from_J = @(J) sqrt(max(1 - 4*real(det(J))/(real(trace(J))^2), 0));

% Safe polarization entropy from coherency
hpol_from_J = @(J) local_hpol_from_J(J);

%% -----------------------------
% Construct synthetic propagated coherency state J_out(t)
% -----------------------------
% The synthetic trajectory is prescribed directly at the coherency level.
% This keeps the demonstration focused on the measurement framework rather
% than on detailed receiver tracking-loop simulation.
Jout = zeros(2,2,N);
S0_true = zeros(N,1);
DoP_true = zeros(N,1);
Hpol_true = zeros(N,1);
S_true = zeros(N,4); % [S0 S1 S2 S3]

% Base total intensity with mild structured variability
S0_series = S0_base * ( ...
    1 ...
    + 0.06*sin(2*pi*t/20) ...
    + 0.025*sin(2*pi*t/5) ...
    + 0.010*movavg(randn(N,1), 21) );
S0_series = max(S0_series, 0.75);

% Smooth correlated perturbations
eta1 = movavg(randn(N,1), 31);
eta2 = movavg(randn(N,1), 31);
eta3 = movavg(randn(N,1), 31);

eta1 = eta1 / max(abs(eta1));
eta2 = eta2 / max(abs(eta2));
eta3 = eta3 / max(abs(eta3));

% -------------------------------------------------------------------------
% Regime I: baseline unitary trajectory
% -------------------------------------------------------------------------
for k = idx1
    th1 = 2*pi*(k-1)/280;
    th2 = 2*pi*(k-1)/190 + 0.4;

    q = [ ...
        0.46*cos(th1) + 0.05*eta1(k); ...
        0.40*sin(th1) + 0.05*cos(th2) + 0.04*eta2(k); ...
        0.30*sin(th2) + 0.08*cos(0.7*th1) + 0.03*eta3(k) ...
        ];

    q = DoP0 * q / max(norm(q), 1e-10);

    Svec = S0_series(k) * q;
    Jk = J_from_stokes(S0_series(k), Svec);
    Jout(:,:,k) = herm(Jk);
    S_true(k,:) = stokes_from_J(Jout(:,:,k));
end

% -------------------------------------------------------------------------
% Regime II: still unitary, but through a more receiver-sensitive region
% -------------------------------------------------------------------------
for k = idx2
    kk = k - idx2(1) + 1;
    th1 = 2*pi*kk/180 + 0.8;
    th2 = 2*pi*kk/130 + 1.1;

    q = [ ...
        0.72*cos(th1) + 0.05*eta1(k); ...
        0.62*sin(th1) + 0.08*cos(th2) + 0.04*eta2(k); ...
        0.12*sin(th2) + 0.05*cos(0.5*th1) + 0.03*eta3(k) ...
        ];

    q = DoP0 * q / max(norm(q), 1e-10);

    Svec = S0_series(k) * q;
    Jk = J_from_stokes(S0_series(k), Svec);
    Jout(:,:,k) = herm(Jk);
    S_true(k,:) = stokes_from_J(Jout(:,:,k));
end

% -------------------------------------------------------------------------
% Regime III: effective mixing / depolarization
% -------------------------------------------------------------------------
% In this regime the state is generated as a convex mixture of two locally
% unitary states, providing a controlled synthetic example of reduced DoP
% and increased polarization entropy.
rampLen = 45;

for k = idx3
    th1 = 2*pi*(k-1)/180;
    th2 = 2*pi*(k-1)/140 + 0.8;

    q1 = [ ...
        0.62*cos(th1) + 0.06*eta1(k); ...
        0.48*sin(th1) + 0.06*cos(th2); ...
        0.32*sin(th2) + 0.18*cos(0.8*th1) ...
        ];

    q2 = [ ...
        0.50*cos(th2 + 0.5) - 0.05*eta2(k); ...
        0.58*sin(th2) + 0.05*cos(1.3*th1); ...
       -0.28*sin(th1) + 0.22*cos(th2) ...
        ];

    q1 = DoP0 * q1 / max(norm(q1), 1e-10);
    q2 = DoP0 * q2 / max(norm(q2), 1e-10);

    J1 = J_from_stokes(S0_series(k), S0_series(k)*q1);
    J2 = J_from_stokes(S0_series(k), S0_series(k)*q2);

    w = 0.55 + 0.22*sin(2*pi*(k-idx3(1))/220) + 0.08*eta1(k);
    w = min(max(w, 0.08), 0.92);

    Jmix = w*J1 + (1-w)*J2;

    if k <= idx3(1) + rampLen - 1
        a = (k - idx3(1)) / max(rampLen-1,1);
        Jprev = Jout(:,:,idx3(1)-1);
        Jk = (1-a)*Jprev + a*Jmix;
    else
        Jk = Jmix;
    end

    Jout(:,:,k) = herm(Jk);
    S_true(k,:) = stokes_from_J(Jout(:,:,k));
end

% Compute invariants of true Jout
for k = 1:N
    Jk = herm(Jout(:,:,k));
    S0_true(k)  = real(trace(Jk));
    DoP_true(k) = dop_from_J(Jk);
    Hpol_true(k)= hpol_from_J(Jk);
end

% Normalized Stokes vector components for panel (a)
p_true = zeros(N,3);
for k = 1:N
    if S_true(k,1) > 0
        p_true(k,:) = S_true(k,2:4) / S_true(k,1);
    else
        p_true(k,:) = [0 0 0];
    end
end

%% -----------------------------
% Define two different receiver response matrices
% -----------------------------
% These two matrices intentionally represent distinct but well-behaved
% receiver mappings so that scalar observables differ while de-embedded
% second-order quantities can still be compared.
PiA = [1.00, 0.08*exp(1i*0.30);
       0.05*exp(-1i*0.20), 0.92*exp(1i*0.10)];

PiB = [0.88*exp(1i*0.15), 0.14*exp(-1i*0.40);
       0.11*exp(1i*0.55), 1.05*exp(-1i*0.05)];

NA = [noisePowA, 0.0015*exp(1i*0.2);
      0.0015*exp(-1i*0.2), noisePowA];

NB = [noisePowB, 0.0020*exp(-1i*0.4);
      0.0020*exp(1i*0.4), noisePowB];

%% -----------------------------
% Generate receiver covariance matrices
% -----------------------------
% A mild common scalar modulation alpha^2 is included to mimic benign
% common-mode effects that rescale covariance without erasing its
% underlying cross-structure.
RyA = zeros(2,2,N);
RyB = zeros(2,2,N);

alpha2 = 1 + 0.05*sin(2*pi*t/9) + 0.01*randn(N,1);
alpha2 = max(alpha2, 0.8);

for k = 1:N
    RyA(:,:,k) = herm(alpha2(k) * PiA * Jout(:,:,k) * PiA' + NA);
    RyB(:,:,k) = herm(alpha2(k) * PiB * Jout(:,:,k) * PiB' + NB);
end

%% -----------------------------
% Scalar receiver observables
% -----------------------------
% The scalar proxy uses a single receiver-channel power observable and is
% included specifically to illustrate receiver dependence.
PA = squeeze(real(RyA(1,1,:)));
PB = squeeze(real(RyB(1,1,:)));

S4A = nan(N,1);
S4B = nan(N,1);

halfWin = floor(winLen/2);
for k = 1:N
    i1 = max(1, k-halfWin);
    i2 = min(N, k+halfWin);
    xa = PA(i1:i2);
    xb = PB(i1:i2);

    ma = mean(xa);
    mb = mean(xb);

    if ma > 0
        S4A(k) = sqrt(max(mean(xa.^2) - ma^2, 0)) / ma;
    end
    if mb > 0
        S4B(k) = sqrt(max(mean(xb.^2) - mb^2, 0)) / mb;
    end
end

S4A_plot = movavg(S4A, smoothWin);
S4B_plot = movavg(S4B, smoothWin);

%% -----------------------------
% De-embedding
% -----------------------------
% This is the core model-based inversion step used to recover a
% receiver-de-embedded estimate of the propagated coherency state.
JhatA = zeros(2,2,N);
JhatB = zeros(2,2,N);

for k = 1:N
    JhatA(:,:,k) = herm(PiA \ (RyA(:,:,k) - NA) / PiA');
    JhatB(:,:,k) = herm(PiB \ (RyB(:,:,k) - NB) / PiB');
end

S0hatA = zeros(N,1); S0hatB = zeros(N,1);
DoPhatA = zeros(N,1); DoPhatB = zeros(N,1);
HhatA   = zeros(N,1); HhatB   = zeros(N,1);

for k = 1:N
    JA = herm(JhatA(:,:,k));
    JB = herm(JhatB(:,:,k));

    S0hatA(k) = real(trace(JA));
    S0hatB(k) = real(trace(JB));

    DoPhatA(k) = dop_from_J(JA);
    DoPhatB(k) = dop_from_J(JB);

    HhatA(k) = hpol_from_J(JA);
    HhatB(k) = hpol_from_J(JB);
end

% Smoothed display quantities
S0_true_plot   = movavg(S0_true, smoothWin);
DoP_true_plot  = movavg(DoP_true, smoothWin);
Hpol_true_plot = movavg(Hpol_true, smoothWin);

p1_plot = movavg(p_true(:,1), smoothWin);
p2_plot = movavg(p_true(:,2), smoothWin);
p3_plot = movavg(p_true(:,3), smoothWin);

S0hatA_plot  = movavg(S0hatA, smoothWin);
S0hatB_plot  = movavg(S0hatB, smoothWin);
DoPhatA_plot = movavg(DoPhatA, smoothWin);
DoPhatB_plot = movavg(DoPhatB, smoothWin);
HhatA_plot   = movavg(HhatA, smoothWin);
HhatB_plot   = movavg(HhatB, smoothWin);

%% -----------------------------
% Useful diagnostics for caption / text
% -----------------------------
scalarMismatch_mean = mean(abs(S4A_plot - S4B_plot), 'omitnan');
recoveryErrA_rms = sqrt(mean((S0hatA - S0_true).^2));
recoveryErrB_rms = sqrt(mean((S0hatB - S0_true).^2));
DoPRecErrA_rms = sqrt(mean((DoPhatA - DoP_true).^2));
DoPRecErrB_rms = sqrt(mean((DoPhatB - DoP_true).^2));
HRecErrA_rms = sqrt(mean((HhatA - Hpol_true).^2));
HRecErrB_rms = sqrt(mean((HhatB - Hpol_true).^2));

%% -----------------------------
% Figure creation
% -----------------------------
figH = figure('Color','w', ...
              'Units','inches', ...
              'Position',[0.4 0.4 16.0 10.6], ...
              'PaperPositionMode','auto', ...
              'InvertHardcopy','off');

tlo = tiledlayout(figH,2,2,'TileSpacing','compact','Padding','compact');

x1 = t(idx1(end));
x2 = t(idx2(end));

% Marker locations used only in panel (d)
mkD_A = 8:45:N;
mkD_B = 28:45:N;

%% -----------------------------
% Panel (a): propagated coherency trajectory in normalized Stokes space
% -----------------------------
ax1 = nexttile(tlo,1);
plot(t, p1_plot, 'LineWidth', lwMain); hold on;
plot(t, p2_plot, 'LineWidth', lwMain);
plot(t, p3_plot, 'LineWidth', lwMain);
xline(x1, '--k', 'LineWidth', lwXline);
xline(x2, '--k', 'LineWidth', lwXline);
xlabel('Time (s)', 'FontSize', fsAx, 'FontName', fontName);
ylabel('Normalized Stokes components', 'FontSize', fsAx, 'FontName', fontName);
title('(a) Synthetic propagated coherency state', 'FontSize', fsTit, ...
    'FontWeight','bold', 'FontName', fontName);
legend('p_1','p_2','p_3', 'Location','southeast', 'FontSize', fsLeg, ...
    'Box','off', 'FontName', fontName);
grid on; box on;
set(ax1, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);
ylim([-0.85 0.85]);

yl = ylim; yr = yl(2)-yl(1);
text(mean(t(idx1)), yl(2)-0.12*yr, 'I: unitary', ...
    'HorizontalAlignment','center', 'FontSize', fsAnn, ...
    'FontWeight','bold', 'FontName', fontName, ...
    'BackgroundColor','w', 'Margin',1.5);
text(mean(t(idx2)), yl(2)-0.12*yr, 'II: proj.-sensitive', ...
    'HorizontalAlignment','center', 'FontSize', fsAnn, ...
    'FontWeight','bold', 'FontName', fontName, ...
    'BackgroundColor','w', 'Margin',1.5);
text(mean(t(idx3)), yl(2)-0.12*yr, 'III: mixing', ...
    'HorizontalAlignment','center', 'FontSize', fsAnn, ...
    'FontWeight','bold', 'FontName', fontName, ...
    'BackgroundColor','w', 'Margin',1.5);

%% -----------------------------
% Panel (b): receiver-dependent scalar observables
% -----------------------------
ax2 = nexttile(tlo,2);
plot(t, S4A_plot, 'LineWidth', lwMain); hold on;
plot(t, S4B_plot, 'LineWidth', lwMain);
xline(x1, '--k', 'LineWidth', lwXline);
xline(x2, '--k', 'LineWidth', lwXline);
xlabel('Time (s)', 'FontSize', fsAx, 'FontName', fontName);
ylabel('Windowed scalar S_4 proxy', 'FontSize', fsAx, 'FontName', fontName);
title('(b) Receiver-dependent scalar observables', 'FontSize', fsTit, ...
    'FontWeight','bold', 'FontName', fontName);
legend('Receiver A','Receiver B', 'Location','northwest', 'FontSize', fsLeg, ...
    'Box','off', 'FontName', fontName);
grid on; box on;
set(ax2, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);

yl = ylim; yr = yl(2)-yl(1);
text(t(48), yl(2)-0.08*yr, sprintf('mean |\\Delta| = %.3f', scalarMismatch_mean), ...
    'FontSize', fsAnn, 'FontName', fontName, ...
    'BackgroundColor','w', 'Margin',1.5);

%% -----------------------------
% Panel (c): de-embedded intensity recovery
% -----------------------------
ax3 = nexttile(tlo,3);
plot(t, S0_true_plot, 'k-', 'LineWidth', lwRef); hold on;
plot(t, S0hatA_plot, '--', 'LineWidth', lwMain, 'Color', [0 0.4470 0.7410]);
plot(t, S0hatB_plot, ':',  'LineWidth', lwMain, 'Color', [0.8500 0.3250 0.0980]);

xline(x1, '--k', 'LineWidth', lwXline);
xline(x2, '--k', 'LineWidth', lwXline);

xlabel('Time (s)', 'FontSize', fsAx, 'FontName', fontName);
ylabel('Intensity / tr(J)', 'FontSize', fsAx, 'FontName', fontName);
title('(c) Receiver-de-embedded intensity recovery', 'FontSize', fsTit, ...
    'FontWeight','bold', 'FontName', fontName);
legend('True tr(J_{out})','Recovered A','Recovered B', ...
    'Location','northwest', 'FontSize', fsLeg, 'Box','off', ...
    'FontName', fontName);
grid on; box on;
set(ax3, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);
ylim([0.8 1.2]);

yl = ylim; yr = yl(2)-yl(1);
text(t(8), yl(1)+0.06*yr, ...
    sprintf('RMSE(A)=%.3f, RMSE(B)=%.3f', recoveryErrA_rms, recoveryErrB_rms), ...
    'FontSize', fsAnn, 'FontName', fontName, ...
    'BackgroundColor','w', 'Margin',1.5);

%% -----------------------------
% Panel (d): invariant diagnostics, true versus recovered
% -----------------------------
ax4 = nexttile(tlo,4);

yyaxis left
hDoPtrue = plot(t, DoP_true_plot, '-', ...
    'Color', [0 0 0], 'LineWidth', lwRef); hold on;

hDoPA = plot(t(mkD_A), DoPhatA_plot(mkD_A), 'o', ...
    'MarkerSize', mkMed, ...
    'LineStyle', 'none', ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerEdgeColor', [0 0.4470 0.7410]);

hDoPB = plot(t(mkD_B), DoPhatB_plot(mkD_B), 's', ...
    'MarkerSize', mkSmall, ...
    'LineStyle', 'none', ...
    'MarkerFaceColor', [0.3010 0.7450 0.9330], ...
    'MarkerEdgeColor', [0.3010 0.7450 0.9330]);

ylabel('DoP', 'FontSize', fsAx, 'Color', [0 0.4470 0.7410], 'FontName', fontName);
ylim([0.1 0.95]);
ax4.YColor = [0 0.4470 0.7410];

yyaxis right
hHtrue = plot(t, Hpol_true_plot, '-', ...
    'Color', [0.35 0.35 0.35], 'LineWidth', lwRef); hold on;

hHA = plot(t(mkD_A), HhatA_plot(mkD_A), '^', ...
    'MarkerSize', mkMed, ...
    'LineStyle', 'none', ...
    'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'MarkerEdgeColor', [0.8500 0.3250 0.0980]);

hHB = plot(t(mkD_B), HhatB_plot(mkD_B), 'd', ...
    'MarkerSize', mkSmall, ...
    'LineStyle', 'none', ...
    'MarkerFaceColor', [0.9290 0.6940 0.1250], ...
    'MarkerEdgeColor', [0.9290 0.6940 0.1250]);

ylabel('H_{pol}', 'FontSize', fsAx, 'Color', [0.8500 0.3250 0.0980], 'FontName', fontName);
ax4.YColor = [0.8500 0.3250 0.0980];
ylim([0.1 0.8]);

xline(x1, '--k', 'LineWidth', lwXline);
xline(x2, '--k', 'LineWidth', lwXline);
xlabel('Time (s)', 'FontSize', fsAx, 'FontName', fontName);
title('(d) Invariant diagnostics: true and recovered', 'FontSize', fsTit, ...
    'FontWeight','bold', 'FontName', fontName);
grid on; box on;
set(ax4, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);

legend([hDoPtrue, hDoPA, hDoPB, hHtrue, hHA, hHB], ...
    {'DoP true','DoP A','DoP B','H true','H A','H B'}, ...
    'Location','southoutside', 'NumColumns', 3, ...
    'FontSize', fsLeg, 'Box','off', 'FontName', fontName);

%% -----------------------------
% Optional console summary
% -----------------------------
fprintf('Synthetic experiment completed.\n');
fprintf('Regime I  : baseline unitary evolution.\n');
fprintf('Regime II : still unitary, but more receiver-sensitive.\n');
fprintf('Regime III: mixing / depolarization.\n');
fprintf('Receiver A and B produce different scalar metrics, but de-embedded tr(J) should agree closely.\n');
fprintf('RMSE tr(J): A = %.5f, B = %.5f\n', recoveryErrA_rms, recoveryErrB_rms);
fprintf('RMSE DoP : A = %.5f, B = %.5f\n', DoPRecErrA_rms, DoPRecErrB_rms);
fprintf('RMSE Hpol: A = %.5f, B = %.5f\n', HRecErrA_rms, HRecErrB_rms);

%% -----------------------------
% Save figure and synthetic data
% -----------------------------
scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    [scriptDir, ~, ~] = fileparts(scriptFullPath);
end

outDir = fullfile(scriptDir, 'output_figure4');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figBase = fullfile(outDir, 'Figure4_synthetic_coherency_demo');

set(figH, 'Color', 'w');
drawnow;
pause(0.3);

savefig(figH, [figBase '.fig']);

save([figBase '_data.mat'], ...
    't', 'dt', 'N', 'winLen', 'smoothWin', ...
    'idx1', 'idx2', 'idx3', ...
    'Jout', 'S_true', 'S0_true', 'DoP_true', 'Hpol_true', 'p_true', ...
    'PiA', 'PiB', 'NA', 'NB', 'alpha2', ...
    'RyA', 'RyB', 'PA', 'PB', 'S4A', 'S4B', 'S4A_plot', 'S4B_plot', ...
    'JhatA', 'JhatB', 'S0hatA', 'S0hatB', 'S0hatA_plot', 'S0hatB_plot', ...
    'DoPhatA', 'DoPhatB', 'DoPhatA_plot', 'DoPhatB_plot', ...
    'HhatA', 'HhatB', 'HhatA_plot', 'HhatB_plot', ...
    'scalarMismatch_mean', 'recoveryErrA_rms', 'recoveryErrB_rms', ...
    'DoPRecErrA_rms', 'DoPRecErrB_rms', 'HRecErrA_rms', 'HRecErrB_rms');

pngRes = 600;

try
    exportgraphics(figH, [figBase '.pdf'], ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');
catch ME
    warning('PDF export failed: %s', ME.message);
end

try
    set(figH, 'Renderer', 'painters');
    drawnow; pause(0.3);
    print(figH, [figBase '.eps'], '-depsc2', '-painters');
catch ME
    warning('EPS export failed: %s', ME.message);
end

try
    set(figH, 'Renderer', 'opengl');
    drawnow; pause(0.3);
    print(figH, [figBase '.png'], '-dpng', sprintf('-r%d', pngRes));
catch ME
    warning('PNG export failed: %s', ME.message);
end

try
    set(figH, 'Renderer', 'opengl');
    drawnow; pause(0.3);
    print(figH, [figBase '.tif'], '-dtiff', sprintf('-r%d', pngRes));
catch ME
    warning('TIFF export failed: %s', ME.message);
end

try
    exportgraphics(figH, [figBase '.svg'], ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');
catch ME
    warning('SVG export not available or failed: %s', ME.message);
end

fprintf('\nSaved outputs to:\n%s\n', outDir);
fprintf('Saved files:\n');
fprintf('  %s.fig\n', figBase);
fprintf('  %s.pdf\n', figBase);
fprintf('  %s.png\n', figBase);
fprintf('  %s.eps\n', figBase);
fprintf('  %s.svg (if supported)\n', figBase);
fprintf('  %s.tif (if supported)\n', figBase);
fprintf('  %s_data.mat\n\n', figBase);

%% -----------------------------
% Local function
% -----------------------------
function H = local_hpol_from_J(J)
    % Compute polarization entropy from a Hermitian coherency matrix.
    ev = real(eig(0.5*(J + J')));
    ev = max(ev, 0);
    s = sum(ev);
    if s <= 0
        H = 0;
        return;
    end
    p = ev / s;
    p = p(p > 0);
    H = -sum(p .* log(p));
end