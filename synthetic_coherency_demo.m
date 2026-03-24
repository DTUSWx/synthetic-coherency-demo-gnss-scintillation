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
% 2) Receiver de-embedding can recover receiver-normalized coherency
%    quantities that agree across receivers.
% 3) Basis-invariant quantities such as degree of polarization (DoP) and
%    polarization entropy remain constant under unitary polarization
%    rotation, but change under effective mixing / depolarization.
%
% Output
% ------
% A four-panel figure suitable as a draft Figure 4, plus publication-quality
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
% Author: Tibor Durgonics / ChatGPT draft support
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
idx1 = 1:400;              % Regime I: unitary rotation
idx2 = 401:800;            % Regime II: same field viewed by different receivers
idx3 = 801:1200;           % Regime III: mixing / depolarization

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
lwMain = 1.9;
lwAux  = 1.4;
lwRef  = 2.2;
fsAx   = 12;
fsTit  = 13;
fsLeg  = 10;

%% -----------------------------
% Pauli matrices and helpers
% -----------------------------
I2 = eye(2);
sigma1 = [0 1; 1 0];
sigma2 = [0 -1i; 1i 0];
sigma3 = [1 0; 0 -1];

% Anonymous helpers
herm = @(X) 0.5 * (X + X');
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

% Smooth correlated perturbations to avoid overly deterministic motion
eta1 = movavg(randn(N,1), 31);
eta2 = movavg(randn(N,1), 31);
eta3 = movavg(randn(N,1), 31);

eta1 = eta1 / max(abs(eta1));
eta2 = eta2 / max(abs(eta2));
eta3 = eta3 / max(abs(eta3));

% -------------------------------------------------------------------------
% Regimes I and II: unitary polarization evolution with constant DoP
% We construct a smoothly evolving normalized Stokes vector p(t) with
% nontrivial variation in all three components, then enforce |p| = DoP0.
% -------------------------------------------------------------------------
for k = [idx1 idx2]
    th1 = 2*pi*(k-1)/250;
    th2 = 2*pi*(k-1)/170 + 0.6;

    q = [ ...
        0.58*cos(th1) + 0.10*eta1(k); ...
        0.52*sin(th1) + 0.08*cos(th2) + 0.08*eta2(k); ...
        0.38*sin(th2) + 0.10*cos(0.7*th1) + 0.06*eta3(k) ...
        ];

    % Enforce constant polarization purity in the unitary regime
    qnorm = norm(q);
    if qnorm < 1e-10
        q = [DoP0; 0; 0];
    else
        q = DoP0 * q / qnorm;
    end

    Svec = S0_series(k) * q;
    Jk = J_from_stokes(S0_series(k), Svec);
    Jout(:,:,k) = herm(Jk);

    S = stokes_from_J(Jout(:,:,k));
    S_true(k,:) = S;
end

% -------------------------------------------------------------------------
% Regime III: effective mixing / depolarization
% Mix two locally unitary states with distinct polarization states.
% A smooth ramp is added at the boundary to avoid a hard transition.
% -------------------------------------------------------------------------
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

    % Smoothly varying mixture weight with mild structured variability
    w = 0.55 + 0.22*sin(2*pi*(k-idx3(1))/220) + 0.08*eta1(k);
    w = min(max(w, 0.08), 0.92);

    Jmix = w*J1 + (1-w)*J2;

    % Smooth transition into the mixing regime
    if k <= idx3(1) + rampLen - 1
        a = (k - idx3(1)) / max(rampLen-1,1);
        Jprev = Jout(:,:,idx3(1)-1);
        Jk = (1-a)*Jprev + a*Jmix;
    else
        Jk = Jmix;
    end

    Jout(:,:,k) = herm(Jk);

    S = stokes_from_J(Jout(:,:,k));
    S_true(k,:) = S;
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
% Receiver A: mild imbalance and cross-polar leakage
PiA = [1.00, 0.08*exp(1i*0.30);
       0.05*exp(-1i*0.20), 0.92*exp(1i*0.10)];

% Receiver B: stronger mixing / different gain-phase response
PiB = [0.88*exp(1i*0.15), 0.14*exp(-1i*0.40);
       0.11*exp(1i*0.55), 1.05*exp(-1i*0.05)];

% Noise covariance matrices
NA = [noisePowA, 0.0015*exp(1i*0.2);
      0.0015*exp(-1i*0.2), noisePowA];

NB = [noisePowB, 0.0020*exp(-1i*0.4);
      0.0020*exp(1i*0.4), noisePowB];

%% -----------------------------
% Generate receiver covariance matrices
% -----------------------------
RyA = zeros(2,2,N);
RyB = zeros(2,2,N);

% Common scalar weighting alpha(t)^2 to mimic correlator-level nuisance scaling
alpha2 = 1 + 0.05*sin(2*pi*t/9) + 0.01*randn(N,1);
alpha2 = max(alpha2, 0.8);

for k = 1:N
    RyA(:,:,k) = herm(alpha2(k) * PiA * Jout(:,:,k) * PiA' + NA);
    RyB(:,:,k) = herm(alpha2(k) * PiB * Jout(:,:,k) * PiB' + NB);
end

%% -----------------------------
% Scalar receiver observables
% -----------------------------
% Use channel-1 power as a legacy-like single-channel intensity observable
PA = squeeze(real(RyA(1,1,:)));
PB = squeeze(real(RyB(1,1,:)));

% Local S4-like proxy over a sliding window
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

% Smooth only for display
S4A_plot = movavg(S4A, smoothWin);
S4B_plot = movavg(S4B, smoothWin);

%% -----------------------------
% De-embedding
% -----------------------------
JhatA = zeros(2,2,N);
JhatB = zeros(2,2,N);

for k = 1:N
    JhatA(:,:,k) = herm(PiA \ (RyA(:,:,k) - NA) / PiA');
    JhatB(:,:,k) = herm(PiB \ (RyB(:,:,k) - NB) / PiB');
end

% Recovered invariants
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

% Smooth only for display
S0_true_plot = movavg(S0_true, smoothWin);
DoP_true_plot = movavg(DoP_true, smoothWin);
Hpol_true_plot = movavg(Hpol_true, smoothWin);

p1_plot = movavg(p_true(:,1), smoothWin);
p2_plot = movavg(p_true(:,2), smoothWin);
p3_plot = movavg(p_true(:,3), smoothWin);

S0hatA_plot = movavg(S0hatA, smoothWin);
S0hatB_plot = movavg(S0hatB, smoothWin);
DoPhatA_plot = movavg(DoPhatA, smoothWin);
DoPhatB_plot = movavg(DoPhatB, smoothWin);
HhatA_plot   = movavg(HhatA,   smoothWin);
HhatB_plot   = movavg(HhatB,   smoothWin);

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
% Plot Figure 4
% -----------------------------
figH = figure('Units','pixels','Position',[80 80 1420 900], 'Color','w');

x1 = t(idx1(end));
x2 = t(idx2(end));

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% Panel (a): normalized Stokes evolution
% -------------------------------------------------------------------------
nexttile;
plot(t, p1_plot, 'LineWidth', lwMain); hold on;
plot(t, p2_plot, 'LineWidth', lwMain);
plot(t, p3_plot, 'LineWidth', lwMain);
xline(x1, '--k', 'LineWidth', 1.0);
xline(x2, '--k', 'LineWidth', 1.0);
xlabel('Time (s)', 'FontSize', fsAx);
ylabel('Normalized Stokes components', 'FontSize', fsAx);
title('(a) Synthetic propagated coherency state', 'FontSize', fsTit, 'FontWeight','bold');
legend('p_1','p_2','p_3', 'Location','northeast', 'FontSize', fsLeg, 'Box','off');
grid on; box on;
set(gca, 'FontSize', fsAx);
ylim([-0.8 0.9]);

yl = ylim;
yr = yl(2) - yl(1);
text(mean(t(idx1)), yl(2)-0.08*yr, 'Unitary rotation', ...
    'HorizontalAlignment','center', 'FontSize', fsLeg, ...
    'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);
text(mean(t(idx2)), yl(2)-0.08*yr, 'Same field,\newline different receivers', ...
    'HorizontalAlignment','center', 'FontSize', fsLeg, ...
    'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);
text(mean(t(idx3)), yl(2)-0.08*yr, 'Mixing /\newline depolarization', ...
    'HorizontalAlignment','center', 'FontSize', fsLeg, ...
    'FontWeight','bold', 'BackgroundColor','w', 'Margin',1);

% -------------------------------------------------------------------------
% Panel (b): receiver-dependent scalar observables
% -------------------------------------------------------------------------
nexttile;
plot(t, S4A_plot, 'LineWidth', lwMain); hold on;
plot(t, S4B_plot, 'LineWidth', lwMain);
xline(x1, '--k', 'LineWidth', 1.0);
xline(x2, '--k', 'LineWidth', 1.0);
xlabel('Time (s)', 'FontSize', fsAx);
ylabel('Windowed scalar S_4 proxy', 'FontSize', fsAx);
title('(b) Receiver-dependent scalar observables', 'FontSize', fsTit, 'FontWeight','bold');
legend('Receiver A','Receiver B', 'Location','northwest', 'FontSize', fsLeg, 'Box','off');
grid on; box on;
set(gca, 'FontSize', fsAx);

yl = ylim;
yr = yl(2)-yl(1);
text(t(35), yl(2)-0.10*yr, sprintf('mean |\\Delta| = %.3f', scalarMismatch_mean), ...
    'FontSize', fsLeg, 'BackgroundColor','w', 'Margin',1);

% -------------------------------------------------------------------------
% Panel (c): de-embedded receiver-normalized intensity recovery
% -------------------------------------------------------------------------
nexttile;
plot(t, S0_true_plot, 'k-', 'LineWidth', lwRef); hold on;
plot(t, S0hatA_plot, '--', 'LineWidth', lwMain);
plot(t, S0hatB_plot, ':',  'LineWidth', lwRef-0.3);
xline(x1, '--k', 'LineWidth', 1.0);
xline(x2, '--k', 'LineWidth', 1.0);
xlabel('Time (s)', 'FontSize', fsAx);
ylabel('Intensity / tr(J)', 'FontSize', fsAx);
title('(c) De-embedded receiver-normalized recovery', 'FontSize', fsTit, 'FontWeight','bold');
legend('True tr(J_{out})','Recovered A','Recovered B', ...
    'Location','southwest', 'FontSize', fsLeg, 'Box','off');
grid on; box on;
set(gca, 'FontSize', fsAx);

allY = [S0_true_plot; S0hatA_plot; S0hatB_plot];
ymin = min(allY) - 0.03;
ymax = max(allY) + 0.03;
ylim([ymin ymax]);

yl = ylim;
yr = yl(2)-yl(1);
text(t(35), yl(2)-0.10*yr, ...
    sprintf('RMSE(A)=%.3f, RMSE(B)=%.3f', recoveryErrA_rms, recoveryErrB_rms), ...
    'FontSize', fsLeg, 'BackgroundColor','w', 'Margin',1);

% -------------------------------------------------------------------------
% Panel (d): invariant diagnostics, true and recovered
% -------------------------------------------------------------------------
nexttile;

yyaxis left;
plot(t, DoP_true_plot, 'k-', 'LineWidth', lwRef); hold on;
plot(t, DoPhatA_plot, '--', 'LineWidth', lwMain);
plot(t, DoPhatB_plot, ':', 'LineWidth', lwRef-0.3);
ylabel('DoP', 'FontSize', fsAx);
ylim([0.1 0.95]);

yyaxis right;
plot(t, Hpol_true_plot, 'k-', 'LineWidth', lwRef);
plot(t, HhatA_plot, '--', 'LineWidth', lwMain);
plot(t, HhatB_plot, ':', 'LineWidth', lwRef-0.3);
ylabel('H_{pol}', 'FontSize', fsAx);

xline(x1, '--k', 'LineWidth', 1.0);
xline(x2, '--k', 'LineWidth', 1.0);
xlabel('Time (s)', 'FontSize', fsAx);
title('(d) Invariant diagnostics: true and recovered', 'FontSize', fsTit, 'FontWeight','bold');
grid on; box on;
set(gca, 'FontSize', fsAx);

legend({'DoP true','DoP A','DoP B','H true','H A','H B'}, ...
    'Location','southwest', 'FontSize', fsLeg, 'Box','off');

yyaxis left;
yl = ylim;
yr = yl(2)-yl(1);
text(t(35), yl(2)-0.10*yr, ...
    sprintf('DoP RMSE(A)=%.3f, DoP RMSE(B)=%.3f', DoPRecErrA_rms, DoPRecErrB_rms), ...
    'FontSize', fsLeg, 'BackgroundColor','w', 'Margin',1);

%% -----------------------------
% Optional console summary
% -----------------------------
fprintf('Synthetic experiment completed.\n');
fprintf('Regime I-II: unitary rotation, invariants should remain approximately constant.\n');
fprintf('Regime III: mixing/depolarization, DoP should decrease and H_pol should increase.\n');
fprintf('Receiver A and B produce different scalar metrics, but de-embedded tr(J) should agree.\n');
fprintf('RMSE tr(J): A = %.5f, B = %.5f\n', recoveryErrA_rms, recoveryErrB_rms);
fprintf('RMSE DoP : A = %.5f, B = %.5f\n', DoPRecErrA_rms, DoPRecErrB_rms);
fprintf('RMSE Hpol: A = %.5f, B = %.5f\n', HRecErrA_rms, HRecErrB_rms);

%% -----------------------------
% Save figure and synthetic data in publication-quality formats
% -----------------------------
scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
    scriptBaseName = 'synthetic_coherency_demo';
else
    [scriptDir, scriptBaseName, ~] = fileparts(scriptFullPath);
end

outDir = fullfile(scriptDir, 'output_figure4');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

set(figH, 'Color', 'w');

figBase = fullfile(outDir, 'Figure4_synthetic_coherency_demo');

% Save editable MATLAB figure
savefig(figH, [figBase '.fig']);

% Save synthetic data and diagnostics
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

% Preferred export path
try
    exportgraphics(figH, [figBase '.pdf'], 'ContentType', 'vector');
    exportgraphics(figH, [figBase '.png'], 'Resolution', pngRes);
catch ME
    warning('exportgraphics failed: %s\nFalling back to print.', ME.message);
    print(figH, [figBase '.pdf'], '-dpdf', '-painters');
    print(figH, [figBase '.png'], '-dpng', sprintf('-r%d', pngRes));
end

% EPS export
try
    print(figH, [figBase '.eps'], '-depsc2', '-painters');
catch ME
    warning('EPS export failed: %s', ME.message);
end

% SVG export
try
    exportgraphics(figH, [figBase '.svg'], 'ContentType', 'vector');
catch ME
    warning('SVG export not available or failed: %s', ME.message);
end

% TIFF export
try
    print(figH, [figBase '.tif'], '-dtiff', sprintf('-r%d', pngRes));
catch ME
    warning('TIFF export failed: %s', ME.message);
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
