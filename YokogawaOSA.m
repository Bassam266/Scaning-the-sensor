classdef YokogawaOSA < handle
    % YokogawaOSA
    % Robust MEX-based controller for a Yokogawa AQ6370-series OSA.
    %
    % Main corrections compared with the previous version:
    %   1. Uses the requested constructor timeout.
    %   2. Forces internal trigger, SINGLE sweep mode, and Trace A WRITE mode.
    %   3. Clears old data before every sweep.
    %   4. Throws an error when a sweep or trace transfer fails, instead of
    %      silently leaving the previous spectrum in obj.data.
    %   5. Checks that a real wavelength sweep was returned.
    %   6. Safely closes partially constructed or disconnected objects.
    %   7. mexOSASetRen expects a double argument, not int8 (fixed in
    %      SetRemote and Close).
    %
    % The file must be named exactly:
    %   YokogawaOSA.m

    properties (SetAccess = private, GetAccess = public)
        ethernet = struct('ip', '10.48.20.178', ...
                          'port', '10001', ...
                          'user', 'yokogawa', ...
                          'password', '1234');

        buffer = struct('string', [], 'size', []);

        measCond
        measParams = struct('start', [], 'stop', [], 'center', [], ...
                            'span', [], 'res', [], 'sens', [], ...
                            'avg', [], 'nSamples', [], ...
                            'autoSamples', [], 'sweepSpeed', []);

        data = struct('X', [], 'Y', [], 'Ymean', [], 'Ystd', []);
        units = struct('X', [], 'Y', []);

        activeTrace = 'TRA'
        errCode = 0
        errMsg = ''
        timeout = 60000
        isConnected = false
    end

    methods (Access = public)

        function obj = YokogawaOSA(opts)
            arguments
                opts.ip char = '10.48.20.178'
                opts.port char = '10001'
                opts.user char = 'yokogawa'
                opts.password char = '1234'
                opts.timeout (1,1) double {mustBeNumeric, mustBePositive} = 60000
            end

            obj.ethernet.ip = opts.ip;
            obj.ethernet.port = opts.port;
            obj.ethernet.user = opts.user;
            obj.ethernet.password = opts.password;
            obj.timeout = opts.timeout;

            connectionString = [obj.ethernet.ip ',' ...
                                obj.ethernet.port ',' ...
                                obj.ethernet.user ',' ...
                                obj.ethernet.password];

            try
                obj.errCode = mexOSAComStart(4, connectionString);
            catch ME
                obj.isConnected = false;
                error('YokogawaOSA:ConnectionFailed', ...
                    'mexOSAComStart failed: %s', ME.message);
            end

            if obj.errCode ~= 0
                obj.GetErrorMessage();
                obj.isConnected = false;
                error('YokogawaOSA:ConnectionFailed', ...
                    'Could not connect to the OSA: %s', obj.errMsg);
            end

            obj.isConnected = true;
            obj.SetTimeout(opts.timeout);
        end


        function SetASEDefaults(obj, opts)
            % Configure the 1018-1038 nm Fabry-Perot calibration sweep.

            arguments
                obj
                opts.timeout (1,1) double {mustBeNumeric, mustBePositive} = 60000
                opts.remote (1,1) {mustBeNumericOrLogical} = 1
                opts.activetrace char = 'A'
            end

            obj.SetTimeout(opts.timeout);
            obj.GetID();
            obj.SetRemote(opts.remote);

            % Do not wait for an external trigger left over from an earlier
            % experiment.
            obj.Send(':TRIGGER:STATE OFF');

            % Make each initiation a completed single sweep.
            obj.Send(':INITIATE:SMODE SINGLE');

            obj.GetMeasCond();

            % start, stop, center, span, resolution, sensitivity,
            % averages, number of samples, auto samples, sweep speed
            %
            % center/span are computed explicitly (rather than passed as
            % []) so the sweep configuration is self-consistent even if
            % the instrument was left in a stale CENTER/SPAN state (e.g.
            % SPAN = 0) from a previous session.
            startWl_nm = 1018;
            stopWl_nm  = 1038;
            obj.SetMeasCond(startWl_nm, stopWl_nm, ...
                mean([startWl_nm, stopWl_nm]), stopWl_nm - startWl_nm, ...
                0.1, 0, 1, [], 1, 1);

            obj.SetActiveTrace(opts.activetrace);

            % Active trace and trace attribute are separate OSA settings.
            % A trace left in FIX/MAX/MIN mode will not update normally.
            traceName = ['TR' upper(opts.activetrace)];
            obj.Send([':TRACE:STATE:' traceName ' ON']);
            obj.Send([':TRACE:ATTRIBUTE:' traceName ' WRITE']);

            obj.Send('*CLS');

            % Verify the settings before starting a long motor scan.
            obj.VerifySweepConfiguration();

            % Acquire one real spectrum as a preflight test.
            obj.SweepAndRetrieve();
            obj.GetTraceUnits();
        end


        function SetThinFilmsDefault(obj, opts)
            arguments
                obj
                opts.timeout (1,1) double {mustBeNumeric, mustBePositive} = 60000
                opts.remote (1,1) {mustBeNumericOrLogical} = 1
                opts.activetrace char = 'A'
            end

            obj.SetTimeout(opts.timeout);
            obj.GetID();
            obj.SetRemote(opts.remote);
            obj.Send(':TRIGGER:STATE OFF');
            obj.Send(':INITIATE:SMODE SINGLE');

            obj.GetMeasCond();
            startWl_nm = 780;
            stopWl_nm  = 1200;
            obj.SetMeasCond(startWl_nm, stopWl_nm, ...
                mean([startWl_nm, stopWl_nm]), stopWl_nm - startWl_nm, ...
                2, 3, 1, [], 1, 1);

            obj.SetActiveTrace(opts.activetrace);
            traceName = ['TR' upper(opts.activetrace)];
            obj.Send([':TRACE:STATE:' traceName ' ON']);
            obj.Send([':TRACE:ATTRIBUTE:' traceName ' WRITE']);

            obj.Send('*CLS');
            obj.VerifySweepConfiguration();
            obj.SweepAndRetrieve();
            obj.GetTraceUnits();
        end


        function SetTimeout(obj, timeout_ms)
            obj.AssertConnected();

            obj.errCode = mexOSASetTimeout(timeout_ms / 100);
            obj.CheckError('setting communication timeout');

            obj.timeout = timeout_ms;
        end


        function GetMeasCond(obj)
            obj.AssertConnected();

            [obj.errCode, obj.measCond] = mexOSAGetMeasCond();
            obj.CheckError('reading OSA measurement conditions');

            obj.FillReadableStructure();
        end


        function SetMeasCond(obj, start, stop, center, span, res, sens, ...
                avg, nSamples, autoSamples, sweepSpeed)

            obj.AssertConnected();

            if isempty(obj.measCond)
                obj.GetMeasCond();
            end

            setValue = {start; stop; center; span; res; sens; avg; ...
                        nSamples; autoSamples; sweepSpeed};

            % IMPORTANT: [] means "leave this field as currently read from
            % the instrument," not "set this field to empty." Overwriting
            % blindly here previously replaced center/span (or any other
            % field left as []) with empty values, which the instrument
            % interpreted as SPAN = 0 -- collapsing the sweep to a single
            % wavelength at the center. Only overwrite fields that were
            % actually given a value.
            for k = 1:numel(setValue)
                if ~isempty(setValue{k})
                    obj.measCond{k, 2} = setValue{k};
                end
            end

            obj.errCode = mexOSASetMeasCond(obj.measCond);
            obj.CheckError('setting OSA measurement conditions');

            % Read the conditions back from the instrument. This avoids
            % saving requested values when the OSA rejected or adjusted one.
            obj.GetMeasCond();
        end


        function SetActiveTrace(obj, traceLetter)
            obj.AssertConnected();

            traceLetter = upper(traceLetter);
            if numel(traceLetter) ~= 1 || ~ismember(traceLetter, 'ABCDEFG')
                error('YokogawaOSA:InvalidTrace', ...
                    'Trace must be one letter from A to G.');
            end

            % :TRACE:ACTIVE is a set-only command; the instrument sends no
            % reply. Using SendReceive here waited for a response that
            % never arrived and timed out.
            obj.errCode = mexOSASend([':TRACE:ACTIVE TR' traceLetter]);

            obj.CheckError('selecting the active trace');
            obj.activeTrace = ['TR' traceLetter];
        end


        function SweepAndRetrieve(obj)
            % Perform one completed sweep and retrieve only the new trace.

            obj.AssertConnected();

            % Critical: remove the old MATLAB-side spectrum before sweeping.
            % If the sweep fails, stale data can no longer be saved silently.
            obj.data.X = [];
            obj.data.Y = [];
            obj.data.Ymean = [];
            obj.data.Ystd = [];

            % Reassert the essential state before every acquisition. This
            % protects against manual front-panel changes during a scan.
            obj.Send(':TRIGGER:STATE OFF');
            obj.Send(':INITIATE:SMODE SINGLE');
            obj.Send([':TRACE:STATE:' obj.activeTrace ' ON']);
            obj.Send([':TRACE:ATTRIBUTE:' obj.activeTrace ' WRITE']);

            obj.errCode = mexOSASweepSingle();
            obj.CheckError('performing a single OSA sweep');

            [obj.errCode, newX, newY] = mexOSAGetTrace(obj.activeTrace);
            obj.CheckError('retrieving the newly acquired OSA trace');

            newX = newX(:);
            newY = newY(:);

            if isempty(newX) || isempty(newY)
                error('YokogawaOSA:EmptyTrace', ...
                    'The OSA returned an empty wavelength or intensity array.');
            end

            if numel(newX) ~= numel(newY)
                error('YokogawaOSA:InvalidTrace', ...
                    ['The OSA returned %d wavelength values but %d ' ...
                     'intensity values.'], numel(newX), numel(newY));
            end

            if numel(newX) < 10
                error('YokogawaOSA:TooFewPoints', ...
                    'The OSA returned only %d spectral points.', numel(newX));
            end

            if any(~isfinite(newX)) || any(~isfinite(newY))
                error('YokogawaOSA:InvalidTrace', ...
                    'The OSA trace contains NaN or Inf values.');
            end

            wavelengthSpan_m = max(newX) - min(newX);

            % The ASE default span is 20 nm. A vector stuck at one center
            % wavelength has approximately zero wavelength span.
            if wavelengthSpan_m < 1e-9
                error('YokogawaOSA:NoWavelengthSweep', ...
                    ['The returned wavelength axis spans only %.6g nm. ' ...
                     'The OSA did not return a real wavelength sweep.'], ...
                    wavelengthSpan_m * 1e9);
            end

            obj.data.X = newX;
            obj.data.Y = newY;
        end


        function SweepLoopAndRetrieve(obj, nLoops)
            arguments
                obj
                nLoops (1,1) double {mustBeInteger, mustBePositive}
            end

            allY = [];
            commonX = [];

            for n = 1:nLoops
                fprintf('\tSweep %d/%d... ', n, nLoops);
                obj.SweepAndRetrieve();

                if n == 1
                    commonX = obj.data.X;
                    allY = nan(numel(commonX), nLoops);
                else
                    if numel(obj.data.X) ~= numel(commonX) || ...
                            any(obj.data.X ~= commonX)
                        error('YokogawaOSA:WavelengthAxisChanged', ...
                            'The wavelength axis changed between repeated sweeps.');
                    end
                end

                allY(:, n) = obj.data.Y;
                fprintf('Done.\n');
            end

            obj.data.X = commonX;
            obj.data.Y = allY;
            obj.data.Ymean = mean(allY, 2);
            obj.data.Ystd = std(allY, 0, 2);
        end


        function GetTraceUnits(obj)
            obj.AssertConnected();

            [obj.errCode, Xunit, obj.buffer.size] = ...
                mexOSASendReceive(':UNIT:X?');
            obj.CheckError('reading the OSA X-axis unit');

            [obj.errCode, Yunit, obj.buffer.size] = ...
                mexOSASendReceive(':DISPLAY:TRACE:Y1:UNIT?');
            obj.CheckError('reading the OSA Y-axis unit');

            [obj.units.X, obj.units.Y] = ...
                obj.ConvertXYUnits(Xunit(1), Yunit(1));
        end


        function hfig = PlotSpectrum(obj, opts)
            arguments
                obj
                opts.fignum (1,1) double ...
                    {mustBeInteger, mustBePositive, mustBeNonzero} = 1
                opts.persistenceFlag (1,1) logical = false
            end

            if isempty(obj.data.X) || isempty(obj.data.Y)
                error('YokogawaOSA:NoData', ...
                    'No spectrum is available to plot.');
            end

            if ~opts.persistenceFlag
                hfig = figure(opts.fignum);
                clf(hfig);
            else
                hfig = gcf;
                hold on;
            end

            if isvector(obj.data.Y)
                plot(obj.data.X * 1e9, obj.data.Y);
            else
                plot(obj.data.X * 1e9, obj.data.Ymean);
            end

            xlabel('Wavelength [nm]');

            if isempty(obj.units.Y)
                ylabel('OSA level');
            else
                ylabel(obj.units.Y);
            end

            title('Yokogawa OSA spectrum');
            grid on;
        end


        function Send(obj, cmdString)
            obj.AssertConnected();

            obj.errCode = mexOSASend(cmdString);
            obj.CheckError(['sending command "' cmdString '"']);
        end


        function Receive(obj)
            obj.AssertConnected();

            [obj.errCode, obj.buffer.string, obj.buffer.size] = mexOSAReceive();
            obj.CheckError('receiving data from the OSA');
        end


        function SendReceive(obj, cmdString)
            obj.AssertConnected();

            [obj.errCode, obj.buffer.string, obj.buffer.size] = ...
                mexOSASendReceive(cmdString);

            obj.CheckError(['sending query "' cmdString '"']);
        end


        function GetID(obj)
            obj.SendReceive('*IDN?');
        end


        function GetErrorMessage(obj)
            if obj.errCode ~= 0
                try
                    msg = mexOSAGetLastError();
                    if isempty(msg)
                        % Some failures (e.g. a communication timeout) have
                        % no instrument-side error text. Fall back to the
                        % numeric code so the error is never silently blank.
                        msg = sprintf('(no error text available, errCode=%g)', ...
                            obj.errCode);
                    end
                    obj.errMsg = msg;
                catch
                    obj.errMsg = sprintf( ...
                        'Unknown Yokogawa MEX error (errCode=%g).', obj.errCode);
                end
            else
                obj.errMsg = '';
            end
        end


        function SetRemote(obj, remote)
            arguments
                obj
                remote (1,1) {mustBeNumericOrLogical} = 1
            end

            obj.AssertConnected();

            % mexOSASetRen expects a double argument, not int8 - the MEX
            % gateway throws "Input argument error" on non-double input.
            obj.errCode = mexOSASetRen(double(logical(remote)));
            obj.CheckError('changing OSA remote/local state');
        end


        function VerifySweepConfiguration(obj)
            % Confirm the three settings that prevent stale fixed traces.

            triggerState = obj.QueryNumeric(':TRIGGER:STATE?');
            if triggerState ~= 0
                error('YokogawaOSA:ExternalTriggerEnabled', ...
                    'External trigger is still enabled (state %g).', triggerState);
            end

            sweepMode = obj.QueryNumeric(':INITIATE:SMODE?');
            if sweepMode ~= 1
                error('YokogawaOSA:WrongSweepMode', ...
                    'OSA sweep mode is %g, but SINGLE mode should return 1.', ...
                    sweepMode);
            end

            traceAttribute = ...
                obj.QueryNumeric([':TRACE:ATTRIBUTE:' obj.activeTrace '?']);

            if traceAttribute ~= 0
                error('YokogawaOSA:TraceNotWrite', ...
                    ['%s attribute is %g. WRITE mode must return 0; ' ...
                     'FIX/MAX/MIN traces can return old data.'], ...
                    obj.activeTrace, traceAttribute);
            end

            % Confirm the instrument is actually reporting the start/stop
            % wavelengths just requested, not a leftover CENTER/SPAN
            % configuration from a previous session.
            startWl_nm = obj.QueryNumeric(':SENSE:WAVELENGTH:START?') * 1e9;
            stopWl_nm  = obj.QueryNumeric(':SENSE:WAVELENGTH:STOP?') * 1e9;

            if (stopWl_nm - startWl_nm) < 1
                error('YokogawaOSA:StaleRangeConfiguration', ...
                    ['Instrument reports a %.4f-%.4f nm range after ' ...
                     'configuration. This usually means it is still in ' ...
                     'CENTER/SPAN mode from a previous session rather ' ...
                     'than the requested START/STOP range.'], ...
                    startWl_nm, stopWl_nm);
            end
        end


        function Close(obj)
            if ~obj.isConnected
                return;
            end

            % Close should be safe even during error handling or destructor
            % execution. Do not throw a new exception from here.
            try
                mexOSASetRen(double(0));
            catch
            end

            try
                obj.errCode = mexOSAComEnd();
                if obj.errCode ~= 0
                    obj.GetErrorMessage();
                end
            catch
            end

            obj.isConnected = false;
        end


        function delete(obj)
            try
                obj.Close();
            catch
            end
        end

    end


    methods (Access = private)

        function FillReadableStructure(obj)
            if isempty(obj.measCond) || size(obj.measCond, 1) < 10
                return;
            end

            obj.measParams.start = obj.measCond{1, 2};
            obj.measParams.stop = obj.measCond{2, 2};
            obj.measParams.center = obj.measCond{3, 2};
            obj.measParams.span = obj.measCond{4, 2};
            obj.measParams.res = obj.measCond{5, 2};
            obj.measParams.sens = obj.measCond{6, 2};
            obj.measParams.avg = obj.measCond{7, 2};
            obj.measParams.nSamples = obj.measCond{8, 2};
            obj.measParams.autoSamples = obj.measCond{9, 2};
            obj.measParams.sweepSpeed = obj.measCond{10, 2};
        end


        function AssertConnected(obj)
            if ~obj.isConnected
                error('YokogawaOSA:NotConnected', ...
                    'The Yokogawa OSA is not connected.');
            end
        end


        function CheckError(obj, context)
            if obj.errCode == 0
                obj.errMsg = '';
                return;
            end

            obj.GetErrorMessage();
            error('YokogawaOSA:CommunicationError', ...
                'Error while %s: %s', context, obj.errMsg);
        end


        function value = QueryNumeric(obj, cmdString)
            obj.SendReceive(cmdString);

            response = strtrim(char(obj.buffer.string));
            value = str2double(response);

            if ~isfinite(value)
                error('YokogawaOSA:InvalidQueryResponse', ...
                    'Query %s returned "%s", not a number.', ...
                    cmdString, response);
            end
        end

    end


    methods (Static)

        function [unitXstr, unitYstr] = ConvertXYUnits(unitXcode, unitYcode)

            switch str2double(unitXcode)
                case 0
                    unitXstr = 'Wavelength [m]';
                case 1
                    unitXstr = 'Frequency [Hz]';
                case 2
                    unitXstr = 'Wavenumber [m^{-1}]';
                otherwise
                    unitXstr = 'Unknown X unit';
            end

            switch str2double(unitYcode)
                case 0
                    unitYstr = 'Power [dBm]';
                case 1
                    unitYstr = 'Power [W]';
                case 2
                    unitYstr = 'Power spectral density [dBm/nm]';
                case 3
                    unitYstr = 'Power spectral density [W/nm]';
                otherwise
                    unitYstr = 'Unknown Y unit';
            end
        end

    end
end
