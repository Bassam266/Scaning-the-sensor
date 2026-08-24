%% FP_scan_acquire.m
% Motor + OSA raster scan -- ACQUISITION ONLY. Organized into 6 blocks:
%   1. Settings & scan grid
%   2. OSA connect & configure
%   3. Motor connect & configure
%   4. Integrated raster scan (motor + OSA acquisition)
%   5. Save data
%   6. Cleanup
%
% Output: a .mat file containing lambda, spectra, scan, osaCfg, motorCfg.
% Run FP_calib_analyze.m afterward (pointing it at that file) to analyze
% and visualize the result -- that part has no hardware dependency and can
% be re-run anytime without touching the OSA or motors.
%
% Requires on path: YokogawaOSA.m, KIM101Motor.m, Thorlabs Kinesis DLLs.

clear; clc; close all;

%% ============================== 1. SETTINGS & SCAN GRID ==============================

% --- OSA ---
osaCfg.ip         = '10.48.20.178';
osaCfg.port       = '10001';
osaCfg.user       = 'yokogawa';
osaCfg.password   = '1234';
osaCfg.timeout_ms = 60000;

osaCfg.expectedStart_nm  = 1018;
osaCfg.expectedStop_nm   = 1038;
osaCfg.minSpan_nm        = 10;     % below this, treat as "no real sweep"
osaCfg.maxRetries        = 2;      % re-tries per pixel before giving up
osaCfg.reconnectPause_s  = 3;
osaCfg.identicalLimit    = 3;      % consecutive identical traces -> error

% --- Motor ---
motorCfg.kimSerial   = '97252642';
motorCfg.kinesisPath = 'C:\Program Files\Thorlabs\Kinesis\';
motorCfg.timeout_ms  = 60000;
motorCfg.umPerStep   = 0.020;      % nominal MPIA10 calibration -- VERIFY with
                                    % MotorStepCalibration.m before trusting a
                                    % fine (<10 um) scan step
motorCfg.settle_s    = 0.2;        % pause after every jog

% --- Scan grid ---
% This script does NOT auto-navigate to a corner. Before running it,
% manually jog the stage (Kinesis GUI) to wherever you want the scan to
% START -- that position becomes (0,0). The scan then moves only in the
% +X/+Y direction from there, covering scanFOV_um x scanFOV_um.
scanFOV_um  = 500;   % total field of view
scanStep_um = 10;      % step between pixels

scan.stepUm = scanStep_um;
scan.nx = round(scanFOV_um / scanStep_um) + 1;
scan.ny = scan.nx;
scan.Npts = scan.nx * scan.ny;
scan.x_um = (0:scan.nx - 1) * scanStep_um;   % relative to wherever you started
scan.y_um = (0:scan.ny - 1) * scanStep_um;
scan.jogSteps = round(scanStep_um / motorCfg.umPerStep);

fprintf('Scan grid: %d x %d = %d points, jog size = %d raw steps\n', ...
    scan.nx, scan.ny, scan.Npts, scan.jogSteps);
fprintf('Set Jog Size Fwd/Rev = %d steps in Kinesis, then disconnect Kinesis GUI.\n', ...
    scan.jogSteps);
fprintf(['Make sure the stage is already positioned at your desired scan-start ' ...
    'corner -- this script will NOT move to a corner automatically.\n\n']);

if ~strcmpi(input('Continue with scan? Type y and press Enter: ', 's'), 'y')
    disp('Scan cancelled.');
    return;
end

%% ============================== 2. OSA CONNECT & CONFIGURE ==============================

fprintf('\nConnecting OSA...\n');

osa = YokogawaOSA('ip', osaCfg.ip, 'port', osaCfg.port, ...
    'user', osaCfg.user, 'password', osaCfg.password, ...
    'timeout', osaCfg.timeout_ms);

osa.SetASEDefaults('timeout', osaCfg.timeout_ms);   % 1018-1038 nm, SINGLE sweep, Trace A WRITE

osa.GetID();
fprintf('Connected: %s\n', strtrim(char(osa.buffer.string)));
fprintf('Preflight sweep: %.2f to %.2f nm, %d points\n\n', ...
    min(osa.data.X) * 1e9, max(osa.data.X) * 1e9, numel(osa.data.X));

%% ============================== 3. MOTOR CONNECT & CONFIGURE ==============================

fprintf('Connecting KIM101 motor...\n');

kim = KIM101Motor('serialNumber', motorCfg.kimSerial, 'kinesisPath', motorCfg.kinesisPath, ...
    'timeout_ms', motorCfg.timeout_ms, 'settle_s', motorCfg.settle_s, ...
    'umPerStep', motorCfg.umPerStep);

% STEP 1: set the current physical position as (0,0).
kim.ZeroPosition();

fprintf('Current position set as zero. Scanning outward from here.\n\n');

%% ============================== 4. INTEGRATED RASTER SCAN ==============================
% STEP 2: raster scan outward from the zeroed position (serpentine: row
% direction alternates to save moves). Acquire a verified spectrum at
% every pixel, then jog to the next one.

lambda  = [];
spectra = [];             % (kx, ky, wavelength) once lambda is known
previousY = [];
identicalCount = 0;
scanTimer = tic;

try
    for ky = 1:scan.ny

        if mod(ky, 2) == 1
            xOrder = 1:scan.nx; xDir = +1;
        else
            xOrder = scan.nx:-1:1; xDir = -1;
        end

        for i = 1:numel(xOrder)
            kx = xOrder(i);

            [osa, currentLambda, currentY, previousY, identicalCount] = ...
                AcquireVerifiedSpectrum(osa, osaCfg, previousY, identicalCount);

            if isempty(spectra)
                lambda = currentLambda;
                spectra = nan(scan.nx, scan.ny, numel(lambda));
            elseif numel(currentLambda) ~= numel(lambda) || ...
                    max(abs(currentLambda - lambda)) > 1e-13
                error('FPScan:WavelengthAxisChanged', ...
                    'OSA wavelength axis changed mid-scan.');
            end

            spectra(kx, ky, :) = currentY;

            fprintf('Pixel %d/%d (kx=%d, ky=%d, X=%.0f um, Y=%.0f um) done. Elapsed %.1f min\n', ...
                (ky - 1) * scan.nx + i, scan.Npts, kx, ky, ...
                scan.x_um(kx), scan.y_um(ky), toc(scanTimer) / 60);

            if i < numel(xOrder)
                kim.Jog('x', xDir);
            end
        end

        if ky < scan.ny
            kim.Jog('y', +1);
        end
    end

catch ME
    warning('Scan interrupted: %s', ME.message);
    kim.ReturnToZero();   % STEP 3, even on error: back to the same position it started at
    SafeDeleteOSA(osa);
    kim.Close();
    rethrow(ME);
end

% STEP 3: scan finished normally -- return to the same position it started at.
fprintf('\nScan complete in %.1f min. Returning to zero...\n', toc(scanTimer) / 60);
kim.ReturnToZero();

%% ============================== 5. SAVE DATA ==============================

folder = fullfile(pwd, datestr(now, 'yyyy-mm-dd'));
if ~exist(folder, 'dir'), mkdir(folder); end

fileName = ['calibscan_sensor162_Fov500um_SS10um' datestr(now, 'HHMMSS') '.mat'];
saveFile = fullfile(folder, fileName);

save(saveFile, 'lambda', 'spectra', 'scan', 'osaCfg', 'motorCfg', '-v7.3');
fprintf('Saved:\n%s\n', saveFile);
fprintf('Run FP_calib_analyze.m and point it at this file to analyze the results.\n');

%% ============================== 6. CLEANUP ==============================

SafeDeleteOSA(osa);
kim.Close();
fprintf('Done.\n');


%% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function [osa, lambda, Y, previousY, identicalCount] = ...
        AcquireVerifiedSpectrum(osa, osaCfg, previousY, identicalCount)
    % Sweeps, validates, and retries/reconnects on failure.

    lastError = [];

    for attempt = 1:(osaCfg.maxRetries + 1)
        try
            if isempty(osa) || ~isvalid(osa) || ~osa.isConnected
                osa = YokogawaOSA('ip', osaCfg.ip, 'port', osaCfg.port, ...
                    'user', osaCfg.user, 'password', osaCfg.password, ...
                    'timeout', osaCfg.timeout_ms);
                osa.SetASEDefaults('timeout', osaCfg.timeout_ms);
            end

            osa.SweepAndRetrieve();
            lambda = osa.data.X(:);
            Y = osa.data.Y(:);

            ValidateSpectrum(lambda, Y, osaCfg.expectedStart_nm, ...
                osaCfg.expectedStop_nm, osaCfg.minSpan_nm);

            if ~isempty(previousY) && numel(previousY) == numel(Y) && ...
                    isequaln(previousY, Y)
                identicalCount = identicalCount + 1;
            else
                identicalCount = 0;
            end

            if identicalCount >= osaCfg.identicalLimit
                error('FPScan:RepeatedIdenticalTrace', ...
                    'Same spectrum returned %d times in a row (stale trace).', ...
                    identicalCount + 1);
            end

            previousY = Y;
            return;

        catch ME
            lastError = ME;
            fprintf(2, 'OSA acquisition attempt %d/%d failed: %s\n', ...
                attempt, osaCfg.maxRetries + 1, ME.message);

            if attempt <= osaCfg.maxRetries
                SafeDeleteOSA(osa);
                pause(osaCfg.reconnectPause_s);
                osa = YokogawaOSA('ip', osaCfg.ip, 'port', osaCfg.port, ...
                    'user', osaCfg.user, 'password', osaCfg.password, ...
                    'timeout', osaCfg.timeout_ms);
                osa.SetASEDefaults('timeout', osaCfg.timeout_ms);
            end
        end
    end

    rethrow(lastError);
end

function ValidateSpectrum(lambda, Y, expectedStart_nm, expectedStop_nm, minSpan_nm)
    if isempty(lambda) || isempty(Y)
        error('FPScan:EmptySpectrum', 'The OSA returned an empty spectrum.');
    end
    if numel(lambda) ~= numel(Y)
        error('FPScan:InvalidSpectrumSize', 'Wavelength/intensity length mismatch.');
    end
    if any(~isfinite(lambda)) || any(~isfinite(Y))
        error('FPScan:NonFiniteSpectrum', 'Spectrum contains NaN or Inf.');
    end

    lambda_nm = lambda * 1e9;
    span_nm = max(lambda_nm) - min(lambda_nm);
    if span_nm < minSpan_nm
        error('FPScan:NoWavelengthSweep', ...
            'Wavelength span only %.4g nm; expected ~%.1f nm.', ...
            span_nm, expectedStop_nm - expectedStart_nm);
    end

    tol_nm = 1.0;
    if min(lambda_nm) > expectedStart_nm + tol_nm || max(lambda_nm) < expectedStop_nm - tol_nm
        error('FPScan:WrongRange', 'OSA range %.2f-%.2f nm does not match expected %.2f-%.2f nm.', ...
            min(lambda_nm), max(lambda_nm), expectedStart_nm, expectedStop_nm);
    end
end

function SafeDeleteOSA(osa)
    if isempty(osa), return; end
    try
        if isvalid(osa), delete(osa); end
    catch
    end
end
