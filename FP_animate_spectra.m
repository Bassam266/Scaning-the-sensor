%% FP_animate_spectra.m
% Animate the reflectivity spectrum at every measurement point of a saved
% FP calibration scan alongside the dip-wavelength map, and save the
% animation as a GIF.
%
%   1. Load the saved scan file
%   2. Compute the dip (minimum) wavelength map, same as FP_calib_analyze.m
%   3. Loop over each measurement point (kx, ky): show the map with the
%      current point marked, next to that point's spectrum
%   4. Save each plotted frame into a GIF
%
% Requires on path: none beyond base MATLAB.

clear; clc; close all;

%% ============================== 1. LOAD SAVED SCAN ==============================

[fname, fpath] = uigetfile('calibscan*.mat', 'Select a calibration scan file');
if isequal(fname, 0)
    error('No file selected.');
end
FILE_TO_ANALYZE = fullfile(fpath, fname);

loaded = load(FILE_TO_ANALYZE, 'lambda', 'spectra', 'scan');
lambda = loaded.lambda;
spectra = loaded.spectra;
scan = loaded.scan;

lambda_nm = lambda(:) * 1e9;

fprintf('Loaded: %s\n', FILE_TO_ANALYZE);
fprintf('spectra size: %s, lambda size: %s\n', ...
    mat2str(size(spectra)), mat2str(size(lambda)));

%% ============================== 2. DIP MAP ==============================

dipVal_map = nan(scan.ny, scan.nx);
dipWl_map  = nan(scan.ny, scan.nx);
dipIdx_map = nan(scan.ny, scan.nx);

for ky = 1:scan.ny
    for kx = 1:scan.nx
        pixSpec = squeeze(spectra(kx, ky, :));
        [dp, dpIdx] = min(pixSpec);

        dipVal_map(ky, kx)  = dp;
        dipWl_map(ky, kx)   = lambda_nm(dpIdx);
        dipIdx_map(ky, kx)  = dpIdx;
    end
end

%% ============================== 3. SETTINGS ==============================

GIF_FILE   = fullfile(fpath, 'spectra_animation.gif');
FRAME_TIME = 0.1;   % seconds per frame

yMin = min(spectra(:), [], 'omitnan');
yMax = max(spectra(:), [], 'omitnan');

%% ============================== 4. PLOT MAP + SPECTRUM & BUILD GIF ==============================

fig = figure('Name', 'Spectra animation', 'Position', [100 100 1000 450]);

isFirstFrame = true;
pt = 0;

for ky = 1:scan.ny
    for kx = 1:scan.nx
        pt = pt + 1;

        pixSpec = squeeze(spectra(kx, ky, :));
        dpIdx = dipIdx_map(ky, kx);

        % --- left panel: dip-wavelength map with current point marked ---
        subplot(1, 2, 1);
        imagesc(scan.x_um, scan.y_um, dipWl_map);
        set(gca, 'YDir', 'normal'); axis image; colormap(gca, 'turbo');
        cb = colorbar;
        cb.FontSize = 14;
        cb.Label.String = 'Dip wavelength [nm]';
        cb.Label.FontSize = 14;
        xlabel('x [um]'); ylabel('y [um]'); title('Dip wavelength map [nm]');
        hold on;
        plot(scan.x_um(kx), scan.y_um(ky), 'wo', 'MarkerSize', 10, 'LineWidth', 2);
        hold off;

        % --- right panel: spectrum at current point, with dip marked ---
        subplot(1, 2, 2);
        plot(lambda_nm, pixSpec, 'LineWidth', 1.5); hold on;
        plot(lambda_nm(dpIdx), dipVal_map(ky, kx), 'bv', 'MarkerFaceColor', 'b');
        hold off;
        xlabel('Wavelength [nm]');
        ylabel('Reflectivity');
        ylim([yMin, yMax]);
        grid on;
        legend('Spectrum', 'Dip', 'Location', 'best');
        title(sprintf('Point %d/%d  (kx=%d, ky=%d)', pt, scan.Npts, kx, ky));

        drawnow;

        % --- capture frame and append to GIF ---
        frame = getframe(fig);
        [imind, cm] = rgb2ind(frame2im(frame), 256);

        if isFirstFrame
            imwrite(imind, cm, GIF_FILE, 'gif', 'Loopcount', inf, 'DelayTime', FRAME_TIME);
            isFirstFrame = false;
        else
            imwrite(imind, cm, GIF_FILE, 'gif', 'WriteMode', 'append', 'DelayTime', FRAME_TIME);
        end
    end
end

close(fig);

fprintf('\nGIF saved to: %s\n', GIF_FILE);
