%% synthetic_coherency_sensitivity_demo.m
% Synthetic sensitivity study for the coherency-based GNSS scintillation paper.
%
% Publication-oriented version
% ----------------------------
% This script generates a compact synthetic sensitivity figure designed to
% support the practical claims of the coherency-based measurement model.
%
% Main features
% -------------
% 1) Finite-sample covariance estimation:
%    The receiver covariance is not treated as perfectly known. Instead,
%    it is estimated from a finite number of synthetic complex samples.
%    This introduces realistic estimation uncertainty.
%
% 2) Regularized de-embedding:
%    Receiver de-embedding is performed with a regularized pseudo-inverse
%    to avoid pathological blow-up in strongly ill-conditioned cases.
%
% 3) Hermitian PSD projection:
%    The recovered coherency matrix is projected onto the Hermitian
%    positive-semidefinite cone, ensuring physically admissible estimates.
%
% 4) Comparable normalized error metrics:
%    Recovery errors for intensity, DoP, and polarization entropy are all
%    reported on roughly comparable percent scales so the line panels are
%    easy to interpret.
%
% 5) Structured calibration perturbations:
%    The calibration-mismatch panel uses gain, phase, and cross-coupling
%    perturbations rather than purely iid matrix noise. This makes the
%    calibration sensitivity trend more informative.
%
% 6) Reviewer-friendly heat-map display:
%    The underlying intensity-recovery error is kept in linear percent
%    units, but panel (d) is displayed using log10(1+error) color scaling
%    so the internal gradient remains visible even when the dynamic range
%    is large.
%
% Figure intent
% -------------
% This figure is designed to support three practical points:
%
%   (i)  observability alone is not sufficient for stable recovery;
%   (ii) recovery degrades as the receiver response becomes more
%        ill-conditioned;
%   (iii) recovery also degrades under calibration mismatch and lower
%         effective SNR.
%
% The figure is qualitative rather than universal. It is a compact
% synthetic sensitivity study, not an operational performance bound.
%
% Output basename
% ---------------
%   Figure_5_sensitivity_demo
%
% -------------------------------------------------------------------------

clear; close all; clc;

%% ========================================================================
%  User settings
% ========================================================================

% Number of synthetic coherency epochs in the truth set
N  = 320;
dt = 0.1;
t  = (0:N-1)' * dt;

% Truth-set intensity and purity scale
S0_base = 1.0;
DoP0    = 0.82;

% Finite-sample covariance support per epoch
% ------------------------------------------
% Larger Nsamp reduces Monte Carlo roughness and weakens purely
% sample-driven blow-up, while still allowing conditioning and SNR
% sensitivity to remain visible.
Nsamp = 420;

% Monte Carlo repetitions
nMC_cond = 24;
nMC_cal  = 24;
nMC_snr  = 24;
nMC_heat = 12;

% Random seed for reproducibility
rng(7);

% Sweep values
condVals   = linspace(1.5, 25.0, 11);   % receiver conditioning sweep
calErrVals = linspace(0.0, 0.12, 11);   % 0 to 12 percent calibration perturbation
snrVals_dB = linspace(35, 5, 11);       % effective covariance SNR sweep

% Heat-map grid
condGrid   = linspace(1.5, 25.0, 11);
snrGrid_dB = linspace(35, 5, 11);

% Regularization floor relative to the largest singular value of Pi
% -----------------------------------------------------------------
% This imposes a lower bound on inverted singular values during
% de-embedding. It is mild, but it prevents the most extreme numerical
% amplification in the worst-conditioned cases.
regFloor = 3e-2;

% Heat-map display controls
% -------------------------
% The actual error metric remains linear percent error.
% These settings affect only the visualization of panel (d).
heatContourLevels_pct = [5 10 20 40 60 80 100 120];
heatTickVals_pct      = [0 5 10 20 40 80 120];

% Mild moderation of DoP and entropy scaling for visual comparability
% -------------------------------------------------------------------
% The intensity error tends to dominate naturally. This factor keeps the
% DoP and Hpol curves visible without changing the qualitative trend.
doEntropyScale = 0.85;

%% ========================================================================
%  Plot styling
% ========================================================================

fontName = 'Helvetica';

lwMain = 2.7;
lwAxis = 1.2;

fsTick = 15;
fsAx   = 17;
fsTit  = 18;
fsLeg  = 13;

colS0   = [0.0000 0.4470 0.7410];
colDoP  = [0.8500 0.3250 0.0980];
colHpol = [0.9290 0.6940 0.1250];

%% ========================================================================
%  Pauli matrices and helpers
% ========================================================================

I2     = eye(2);
sigma1 = [0 1; 1 0];
sigma2 = [0 -1i; 1i 0];
sigma3 = [1 0; 0 -1];

herm = @(X) 0.5 * (X + X');
J_from_stokes = @(S0,S) 0.5 * (S0*I2 + S(1)*sigma1 + S(2)*sigma2 + S(3)*sigma3);

% Maximum polarization entropy for a normalized 2x2 coherency matrix
Hmax = log(2);

%% ========================================================================
%  Build synthetic truth set Jout(k)
% ========================================================================
% The truth set is designed to be smooth and interpretable:
% - most of the interval resembles coherent unitary-like evolution;
% - the final segment introduces mixed-state behavior.

Jout      = zeros(2,2,N);
S0_true   = zeros(N,1);
DoP_true  = zeros(N,1);
Hpol_true = zeros(N,1);

eta1 = smoothdata(randn(N,1),'movmean',15);
eta2 = smoothdata(randn(N,1),'movmean',15);
eta3 = smoothdata(randn(N,1),'movmean',15);

eta1 = eta1 / max(abs(eta1));
eta2 = eta2 / max(abs(eta2));
eta3 = eta3 / max(abs(eta3));

S0_series = S0_base * ( ...
    1 ...
    + 0.05*sin(2*pi*t/18) ...
    + 0.02*sin(2*pi*t/4.5) ...
    + 0.01*eta1 );
S0_series = max(S0_series, 0.78);

for k = 1:N
    th1 = 2*pi*(k-1)/170;
    th2 = 2*pi*(k-1)/115 + 0.4;

    q1 = [ ...
        0.60*cos(th1) + 0.04*eta1(k); ...
        0.42*sin(th1) + 0.05*cos(th2) + 0.03*eta2(k); ...
        0.28*sin(th2) + 0.07*cos(0.6*th1) + 0.03*eta3(k)];

    q1 = DoP0 * q1 / max(norm(q1),1e-10);

    % Final segment introduces mixed-state behavior so the truth set is
    % not unrealistically idealized.
    if k > round(0.72*N)
        q2 = [ ...
            0.45*cos(th2 + 0.4); ...
            0.55*sin(th2); ...
           -0.20*sin(th1) + 0.18*cos(th2)];
        q2 = DoP0 * q2 / max(norm(q2),1e-10);

        J1 = J_from_stokes(S0_series(k), S0_series(k)*q1);
        J2 = J_from_stokes(S0_series(k), S0_series(k)*q2);
        w  = 0.65 + 0.15*sin(2*pi*(k-round(0.72*N))/110);
        Jk = w*J1 + (1-w)*J2;
    else
        Jk = J_from_stokes(S0_series(k), S0_series(k)*q1);
    end

    Jk = project_to_hpsd(Jk);
    Jout(:,:,k) = Jk;

    S0_true(k)   = real(trace(Jk));
    DoP_true(k)  = local_dop_from_J(Jk);
    Hpol_true(k) = local_hpol_from_J(Jk);
end

signalScale = mean(S0_true);

%% ========================================================================
%  Reference receiver model
% ========================================================================

Pi_ref = [1.00,               0.08*exp(1i*0.30); ...
          0.05*exp(-1i*0.18), 0.94*exp(1i*0.08)];

% Fixed singular-vector basis for the conditioning family
% -------------------------------------------------------
% This ensures the conditioning sweep changes the singular-value spread
% while keeping the singular-vector directions fixed.
A0 = randn(2) + 1i*randn(2);
B0 = randn(2) + 1i*randn(2);
[U0,~,~] = svd(A0);
[V0,~,~] = svd(B0);

%% ========================================================================
%  Arrays for summary metrics
% ========================================================================

errS0_cond  = zeros(size(condVals));
errDoP_cond = zeros(size(condVals));
errH_cond   = zeros(size(condVals));

errS0_cal   = zeros(size(calErrVals));
errDoP_cal  = zeros(size(calErrVals));
errH_cal    = zeros(size(calErrVals));

errS0_snr   = zeros(size(snrVals_dB));
errDoP_snr  = zeros(size(snrVals_dB));
errH_snr    = zeros(size(snrVals_dB));

errHeat     = zeros(numel(snrGrid_dB), numel(condGrid));

%% ========================================================================
%  (a) Conditioning sweep
% ========================================================================
% Only conditioning is varied here. SNR is kept fixed at a moderate level.

snrCond_dB   = 20;
noisePowCond = signalScale / (10^(snrCond_dB/10));
N_cond = [noisePowCond, 0.15*noisePowCond*exp(1i*0.2); ...
          0.15*noisePowCond*exp(-1i*0.2), noisePowCond];
N_cond = herm(N_cond);

for i = 1:numel(condVals)
    vals = zeros(nMC_cond,3);
    Pi_true = make_matrix_with_condition_fixed_basis(condVals(i), U0, V0);

    for m = 1:nMC_cond
        vals(m,:) = evaluate_recovery_finite_sample( ...
            Jout, S0_true, DoP_true, Hpol_true, ...
            Pi_true, Pi_true, N_cond, N_cond, Nsamp, regFloor, doEntropyScale);
    end

    errS0_cond  = assign_median(errS0_cond, i, vals(:,1));
    errDoP_cond = assign_median(errDoP_cond, i, vals(:,2));
    errH_cond   = assign_median(errH_cond, i, vals(:,3));
end

%% ========================================================================
%  (b) Calibration-error sweep
% ========================================================================
% Only calibration mismatch is varied here. SNR is again fixed at a
% moderate level. The perturbation is structured:
% - gain perturbation on diagonal terms
% - phase perturbation on diagonal terms
% - stronger cross-coupling perturbation off diagonal

snrCal_dB   = 20;
noisePowCal = signalScale / (10^(snrCal_dB/10));
N_cal_true = [noisePowCal, 0.15*noisePowCal*exp(1i*0.2); ...
              0.15*noisePowCal*exp(-1i*0.2), noisePowCal];
N_cal_true = herm(N_cal_true);

for i = 1:numel(calErrVals)
    epsLevel = calErrVals(i);
    vals = zeros(nMC_cal,3);

    for m = 1:nMC_cal
        Pi_true = Pi_ref;
        Pi_est  = Pi_true;

        gainPert  = epsLevel * 0.6 * randn(2,1);
        phasePert = epsLevel * 2.5 * randn(2,1);
        mixPert   = epsLevel * 1.2 * (randn(2) + 1i*randn(2)) / sqrt(2);

        for jj = 1:2
            Pi_est(jj,jj) = Pi_est(jj,jj) * (1 + gainPert(jj)) * exp(1i*phasePert(jj));
        end

        mixScale = max(abs(Pi_true(:)));
        Pi_est(1,2) = Pi_est(1,2) + mixPert(1,2) * mixScale;
        Pi_est(2,1) = Pi_est(2,1) + mixPert(2,1) * mixScale;

        % Prevent nearly singular estimated models
        if rcond(Pi_est) < 1e-6
            Pi_est = Pi_est + 1e-4*eye(2);
        end

        % Noise covariance estimate is also perturbed, but more mildly
        dN = 0.5 * epsLevel * randn(size(N_cal_true));
        N_est = herm(N_cal_true .* (1 + dN));
        N_est(1,1) = max(real(N_est(1,1)), 1e-6);
        N_est(2,2) = max(real(N_est(2,2)), 1e-6);

        vals(m,:) = evaluate_recovery_finite_sample( ...
            Jout, S0_true, DoP_true, Hpol_true, ...
            Pi_true, Pi_est, N_cal_true, N_est, Nsamp, regFloor, doEntropyScale);
    end

    errS0_cal  = assign_median(errS0_cal, i, vals(:,1));
    errDoP_cal = assign_median(errDoP_cal, i, vals(:,2));
    errH_cal   = assign_median(errH_cal, i, vals(:,3));
end

%% ========================================================================
%  (c) Noise / SNR sweep
% ========================================================================
% Only noise level is varied here. Conditioning is held fixed at a
% representative moderate value.

kappaSNR = 8;
Pi_snr = make_matrix_with_condition_fixed_basis(kappaSNR, U0, V0);

for i = 1:numel(snrVals_dB)
    snrDB = snrVals_dB(i);
    noisePow = signalScale / (10^(snrDB/10));
    N_true = [noisePow, 0.15*noisePow*exp(1i*0.2); ...
              0.15*noisePow*exp(-1i*0.2), noisePow];
    N_true = herm(N_true);

    vals = zeros(nMC_snr,3);
    for m = 1:nMC_snr
        vals(m,:) = evaluate_recovery_finite_sample( ...
            Jout, S0_true, DoP_true, Hpol_true, ...
            Pi_snr, Pi_snr, N_true, N_true, Nsamp, regFloor, doEntropyScale);
    end

    errS0_snr  = assign_median(errS0_snr, i, vals(:,1));
    errDoP_snr = assign_median(errDoP_snr, i, vals(:,2));
    errH_snr   = assign_median(errH_snr, i, vals(:,3));
end

%% ========================================================================
%  (d) Combined conditioning-SNR heat map
% ========================================================================
% This panel shows only the intensity error to keep the heat map compact.

for ic = 1:numel(condGrid)
    Pi_true = make_matrix_with_condition_fixed_basis(condGrid(ic), U0, V0);

    for is = 1:numel(snrGrid_dB)
        snrDB = snrGrid_dB(is);
        noisePow = signalScale / (10^(snrDB/10));
        N_true = [noisePow, 0.15*noisePow*exp(1i*0.2); ...
                  0.15*noisePow*exp(-1i*0.2), noisePow];
        N_true = herm(N_true);

        vals = zeros(nMC_heat,1);
        for m = 1:nMC_heat
            out = evaluate_recovery_finite_sample( ...
                Jout, S0_true, DoP_true, Hpol_true, ...
                Pi_true, Pi_true, N_true, N_true, Nsamp, regFloor, doEntropyScale);
            vals(m) = out(1);
        end

        errHeat(is,ic) = median(vals);
    end
end

%% ========================================================================
%  Figure creation
% ========================================================================

figH = figure('Color','w', ...
              'Units','inches', ...
              'Position',[0.4 0.4 16.2 10.4], ...
              'PaperPositionMode','auto', ...
              'InvertHardcopy','off');

tlo = tiledlayout(figH,2,2,'TileSpacing','compact','Padding','compact');

% Common y-limit for panels (a) and (c)
lineMax_AC = max([errS0_cond(:); errDoP_cond(:); errH_cond(:); ...
                  errS0_snr(:);  errDoP_snr(:);  errH_snr(:)]);
lineYmax_AC = max(20, 1.08 * lineMax_AC);

% Separate y-limit for panel (b)
lineMax_B = max([errS0_cal(:); errDoP_cal(:); errH_cal(:)]);
lineYmax_B = max(15, 1.15 * lineMax_B);

% ------------------------------------------------------------------------
% (a) Conditioning sensitivity
% ------------------------------------------------------------------------
ax1 = nexttile(tlo,1);
plot(condVals, errS0_cond, '-',  'LineWidth', lwMain, 'Color', colS0); hold on;
plot(condVals, errDoP_cond,'--', 'LineWidth', lwMain, 'Color', colDoP);
plot(condVals, errH_cond,  ':',  'LineWidth', lwMain, 'Color', colHpol);
xlabel('\kappa(\Pi)', 'FontSize', fsAx, 'FontName', fontName);
ylabel('Normalized recovery error (%)', 'FontSize', fsAx, 'FontName', fontName);
title('(a) Conditioning sensitivity', 'FontSize', fsTit, 'FontWeight','bold', 'FontName', fontName);
legend('Intensity error: tr(J)', 'Purity error: DoP', 'Entropy error: H_{pol}', ...
    'Location','northwest', 'FontSize', fsLeg, 'Box', 'off', 'FontName', fontName);
grid on; box on;
ylim([0 lineYmax_AC]);
set(ax1, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);

% ------------------------------------------------------------------------
% (b) Calibration-error sensitivity
% ------------------------------------------------------------------------
ax2 = nexttile(tlo,2);
plot(100*calErrVals, errS0_cal, '-',  'LineWidth', lwMain, 'Color', colS0); hold on;
plot(100*calErrVals, errDoP_cal,'--', 'LineWidth', lwMain, 'Color', colDoP);
plot(100*calErrVals, errH_cal,  ':',  'LineWidth', lwMain, 'Color', colHpol);
xlabel('Relative calibration perturbation (%)', 'FontSize', fsAx, 'FontName', fontName);
ylabel('Normalized recovery error (%)', 'FontSize', fsAx, 'FontName', fontName);
title('(b) Calibration-error sensitivity', 'FontSize', fsTit, 'FontWeight','bold', 'FontName', fontName);
legend('Intensity error: tr(J)', 'Purity error: DoP', 'Entropy error: H_{pol}', ...
    'Location','northwest', 'FontSize', fsLeg, 'Box', 'off', 'FontName', fontName);
grid on; box on;
ylim([0 lineYmax_B]);
set(ax2, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);

% ------------------------------------------------------------------------
% (c) Noise / SNR sensitivity
% ------------------------------------------------------------------------
ax3 = nexttile(tlo,3);
plot(snrVals_dB, errS0_snr, '-',  'LineWidth', lwMain, 'Color', colS0); hold on;
plot(snrVals_dB, errDoP_snr,'--', 'LineWidth', lwMain, 'Color', colDoP);
plot(snrVals_dB, errH_snr,  ':',  'LineWidth', lwMain, 'Color', colHpol);
set(ax3, 'XDir', 'reverse');
xlabel('Effective covariance SNR (dB)', 'FontSize', fsAx, 'FontName', fontName);
ylabel('Normalized recovery error (%)', 'FontSize', fsAx, 'FontName', fontName);
title('(c) Noise / SNR sensitivity', 'FontSize', fsTit, 'FontWeight','bold', 'FontName', fontName);
legend('Intensity error: tr(J)', 'Purity error: DoP', 'Entropy error: H_{pol}', ...
    'Location','northwest', 'FontSize', fsLeg, 'Box', 'off', 'FontName', fontName);
grid on; box on;
ylim([0 lineYmax_AC]);
set(ax3, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);

% ------------------------------------------------------------------------
% (d) Combined sensitivity heat map
% ------------------------------------------------------------------------
% The underlying metric remains errHeat in linear percent units.
% Only the color mapping is compressed:
%
%    displayHeat = log10(1 + errHeat)
%
% This preserves internal gradient structure and avoids a saturated block.
ax4 = nexttile(tlo,4);

displayHeat = log10(1 + errHeat);

imagesc(condGrid, snrGrid_dB, displayHeat);
set(ax4, 'YDir', 'reverse');
xlabel('\kappa(\Pi)', 'FontSize', fsAx, 'FontName', fontName);
ylabel('Effective covariance SNR (dB)', 'FontSize', fsAx, 'FontName', fontName);
title('(d) Combined sensitivity: tr(J) error (%)', ...
    'FontSize', fsTit, 'FontWeight','bold', 'FontName', fontName);
set(ax4, 'FontSize', fsTick, 'LineWidth', lwAxis, 'FontName', fontName);
box on;
colormap(ax4, turbo);

% Use full transformed range; log compression itself handles the dynamic range
caxis([min(displayHeat(:)) max(displayHeat(:))]);

% Overlay contours in transformed coordinates, corresponding to percent
% levels in the original metric.
hold(ax4, 'on');
contour(ax4, condGrid, snrGrid_dB, displayHeat, ...
    log10(1 + heatContourLevels_pct), ...
    'LineColor', [0.15 0.15 0.15], 'LineWidth', 0.7);
hold(ax4, 'off');

cb = colorbar;

% Colorbar ticks positioned in transformed units, labeled in percent units
cbTicks = log10(1 + heatTickVals_pct);
cbTicks = cbTicks(cbTicks >= min(displayHeat(:)) & cbTicks <= max(displayHeat(:)));
cb.Ticks = cbTicks;
cb.TickLabels = string(round(10.^cbTicks - 1));

cb.Label.String = 'tr(J) recovery error (%)';
cb.Label.FontSize = fsAx;
cb.Label.FontName = fontName;
cb.FontSize = fsTick;
cb.FontName = fontName;

%% ========================================================================
%  Console summary
% ========================================================================

fprintf('Sensitivity study completed.\n');
fprintf('Finite-sample covariance estimation used with Nsamp = %d samples per epoch.\n', Nsamp);
fprintf('Regularized de-embedding floor = %.3g of largest singular value.\n', regFloor);
fprintf('Panel (d) uses log10(1+error) color compression for display only.\n');

%% ========================================================================
%  Save outputs
% ========================================================================

scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    [scriptDir, ~, ~] = fileparts(scriptFullPath);
end

outDir = fullfile(scriptDir, 'output_figure5');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figBase = fullfile(outDir, 'Figure_5_sensitivity_demo');

savefig(figH, [figBase '.fig']);

save([figBase '_data.mat'], ...
    't', 'dt', 'N', 'Nsamp', 'regFloor', 'doEntropyScale', ...
    'heatContourLevels_pct', 'heatTickVals_pct', ...
    'Jout', 'S0_true', 'DoP_true', 'Hpol_true', ...
    'Pi_ref', ...
    'condVals', 'calErrVals', 'snrVals_dB', ...
    'condGrid', 'snrGrid_dB', ...
    'errS0_cond', 'errDoP_cond', 'errH_cond', ...
    'errS0_cal', 'errDoP_cal', 'errH_cal', ...
    'errS0_snr', 'errDoP_snr', 'errH_snr', ...
    'errHeat', 'displayHeat');

pngRes = 600;

try
    exportgraphics(figH, [figBase '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
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
    exportgraphics(figH, [figBase '.svg'], 'ContentType', 'vector', 'BackgroundColor', 'white');
catch ME
    warning('SVG export not available or failed: %s', ME.message);
end

fprintf('\nSaved outputs to:\n%s\n', outDir);
fprintf('Saved files with basename Figure_5_sensitivity_demo\n');

%% ========================================================================
%  Local helper functions
% ========================================================================

function Pi = make_matrix_with_condition_fixed_basis(kappaTarget, U0, V0)
    % Construct a 2x2 complex matrix with a prescribed condition number
    % while keeping singular-vector directions fixed across the sweep.
    s1 = 1.0;
    s2 = max(1.0 / kappaTarget, 1e-8);
    Pi = U0 * diag([s1 s2]) * V0';
    if abs(det(Pi)) < 1e-12
        Pi = Pi + 1e-8 * eye(2);
    end
end

function out = evaluate_recovery_finite_sample( ...
    Jout, S0_true, DoP_true, Hpol_true, ...
    Pi_true, Pi_est, N_true, N_est, Nsamp, regFloor, doEntropyScale)
% Simulate finite-sample covariance estimation and regularized de-embedding.
%
% Returns three normalized error metrics on percent-like scales:
%   out(1): intensity recovery error in percent of mean true tr(J)
%   out(2): DoP recovery error in percent of the [0,1] DoP range
%   out(3): Hpol recovery error in percent of Hmax for a 2x2 coherency state

    N = size(Jout,3);

    S0hat  = zeros(N,1);
    DoPhat = zeros(N,1);
    Hhat   = zeros(N,1);

    for k = 1:N
        % True receiver covariance
        Ry_true = Pi_true * Jout(:,:,k) * Pi_true' + N_true;
        Ry_true = 0.5 * (Ry_true + Ry_true');

        % Finite-sample covariance estimate
        Y = generate_complex_gaussian_samples(Ry_true, Nsamp);
        Ry_hat = (Y * Y') / Nsamp;
        Ry_hat = 0.5 * (Ry_hat + Ry_hat');

        % De-embedding + PSD projection
        Jhat = deembed_regularized(Ry_hat, Pi_est, N_est, regFloor);
        Jhat = project_to_hpsd(Jhat);

        S0hat(k)  = real(trace(Jhat));
        DoPhat(k) = local_dop_from_J(Jhat);
        Hhat(k)   = local_hpol_from_J(Jhat);
    end

    rmseS0    = sqrt(mean((S0hat - S0_true).^2));
    errS0_pct = 100 * rmseS0 / mean(S0_true);

    rmseDoP    = sqrt(mean((DoPhat - DoP_true).^2));
    errDoP_pct = 100 * doEntropyScale * rmseDoP;

    rmseH    = sqrt(mean((Hhat - Hpol_true).^2));
    errH_pct = 100 * doEntropyScale * rmseH / log(2);

    out = [errS0_pct, errDoP_pct, errH_pct];
end

function Jhat = deembed_regularized(Ry_hat, Pi_est, N_est, regFloor)
    % Regularized pseudo-inverse based de-embedding.
    [U,S,V] = svd(Pi_est, 'econ');
    s = diag(S);
    sReg = max(s, regFloor * max(s));
    Pi_pinv_reg = V * diag(1 ./ sReg) * U';
    Jhat = Pi_pinv_reg * (Ry_hat - N_est) * Pi_pinv_reg';
    Jhat = 0.5 * (Jhat + Jhat');
end

function Y = generate_complex_gaussian_samples(R, Nsamp)
    % Generate circular complex Gaussian samples with covariance R.
    R = 0.5 * (R + R');
    [V,D] = eig(R);
    d = real(diag(D));
    d(d < 0) = 0;
    L = V * diag(sqrt(d));

    Z = (randn(2,Nsamp) + 1i*randn(2,Nsamp)) / sqrt(2);
    Y = L * Z;
end

function Jpsd = project_to_hpsd(J)
    % Project a matrix onto the Hermitian positive-semidefinite cone.
    J = 0.5 * (J + J');
    [V,D] = eig(J);
    d = real(diag(D));
    d(d < 0) = 0;
    Jpsd = V * diag(d) * V';
    Jpsd = 0.5 * (Jpsd + Jpsd');
end

function d = local_dop_from_J(J)
    % Degree of polarization from a 2x2 coherency matrix.
    tj = real(trace(J));
    if tj <= 0
        d = 0;
        return;
    end
    val = 1 - 4 * real(det(J)) / (tj^2);
    d = sqrt(max(val,0));
end

function H = local_hpol_from_J(J)
    % Polarization entropy from a 2x2 coherency matrix.
    ev = real(eig(0.5 * (J + J')));
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

function arr = assign_median(arr, idx, vals)
    % Small helper for readability in sweep loops
    arr(idx) = median(vals);
end