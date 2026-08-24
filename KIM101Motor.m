classdef KIM101Motor < handle
    % KIM101Motor
    % Controller class for a Thorlabs KIM101 + MPIA10/PIA-series piezo
    % inertia actuator, driving two channels (X/Y).
    %
    % VERIFIED against Thorlabs' own example (Motion_Control_Examples,
    % Matlab/KCubes/KIM101/KIM101.m): Connect, WaitForSettingsInitialized,
    % StartPolling, EnableDevice, SetPositionAs, Jog, StopPolling,
    % Disconnect. These are used exactly as documented.
    %
    % BEST-EFFORT (not found in an official example -- wrapped so they
    % fail safely and tell you how to fix them): GetJogParams/SetJogParams
    % (jog step size), GetPosition (controller's own step counter). If
    % these throw on your installed Kinesis version, run
    %   methods(obj.device)
    % to find the correct names and update the two methods flagged below.
    %
    % IMPORTANT -- READ BEFORE TRUSTING ANY POSITION FROM THIS CLASS:
    % This actuator is OPEN-LOOP: no encoder, no limit switches, no real
    % position sensor. Every "position" this class reports (obj.pos, and
    % GetControllerPosition) is a SOFTWARE STEP COUNT, not a physical
    % measurement. Per Thorlabs' own datasheet, actual step size varies
    % ~20-30% run to run "due to component variance, change of direction,
    % and application condition" -- and steps can be lost entirely with
    % no way for software to detect it. Use MotorStepCalibration.m to
    % check umPerStep empirically before trusting fine-step scans, and
    % treat obj.pos as an estimate, not ground truth.

    properties (SetAccess = private, GetAccess = public)
        serialNumber
        kinesisPath
        device
        ch1
        ch2
        jogIncrease
        jogDecrease
        timeout_ms  = 60000
        settle_s    = 0.2
        umPerStep   = 0.020     % nominal um per raw controller step -- CALIBRATE THIS, don't trust it
        isConnected = false

        pos = struct('x', 0, 'y', 0)   % SOFTWARE step counter, relative to last ZeroPosition() call
    end

    methods (Access = public)

        function obj = KIM101Motor(opts)
            arguments
                opts.serialNumber char
                opts.kinesisPath  char   = 'C:\Program Files\Thorlabs\Kinesis\'
                opts.timeout_ms   (1,1) double {mustBePositive} = 60000
                opts.settle_s     (1,1) double {mustBeNonnegative} = 0.2
                opts.umPerStep    (1,1) double {mustBePositive} = 0.020
            end

            obj.serialNumber = opts.serialNumber;
            obj.kinesisPath  = opts.kinesisPath;
            obj.timeout_ms   = opts.timeout_ms;
            obj.settle_s     = opts.settle_s;
            obj.umPerStep    = opts.umPerStep;

            try
                obj.Connect();
            catch ME
                obj.isConnected = false;
                error('KIM101Motor:ConnectionFailed', ...
                    'Could not connect to KIM101 %s: %s', obj.serialNumber, ME.message);
            end
        end


        function Connect(obj)
            NET.addAssembly([obj.kinesisPath 'Thorlabs.MotionControl.DeviceManagerCLI.dll']);
            NET.addAssembly([obj.kinesisPath 'Thorlabs.MotionControl.GenericMotorCLI.dll']);
            motCLI = NET.addAssembly([obj.kinesisPath 'Thorlabs.MotionControl.KCube.InertialMotorCLI.dll']);

            import Thorlabs.MotionControl.DeviceManagerCLI.*
            import Thorlabs.MotionControl.GenericMotorCLI.*
            import Thorlabs.MotionControl.KCube.InertialMotorCLI.*

            DeviceManagerCLI.BuildDeviceList();

            obj.device = KCubeInertialMotor.CreateKCubeInertialMotor(obj.serialNumber);
            obj.device.Connect(obj.serialNumber);
            obj.device.WaitForSettingsInitialized(5000);
            obj.device.StartPolling(250);
            obj.device.EnableDevice();
            pause(1);   % Thorlabs' own example pauses here to let the device finish enabling

            channelsHandle = motCLI.AssemblyHandle.GetType( ...
                'Thorlabs.MotionControl.KCube.InertialMotorCLI.InertialMotorStatus+MotorChannels');
            channelsEnums = channelsHandle.GetEnumValues();

            jogDirectionHandle = motCLI.AssemblyHandle.GetType( ...
                'Thorlabs.MotionControl.KCube.InertialMotorCLI.InertialMotorJogDirection');

            obj.ch1 = channelsEnums.GetValue(0);
            obj.ch2 = channelsEnums.GetValue(1);
            obj.jogIncrease = System.Enum.Parse(jogDirectionHandle, 'Increase');
            obj.jogDecrease = System.Enum.Parse(jogDirectionHandle, 'Decrease');

            obj.isConnected = true;
        end


        function ZeroPosition(obj)
            % Resets the SOFTWARE step counter to zero at the CURRENT
            % physical position. Does not move the stage. Does not
            % establish a real physical reference -- see the class-level
            % note on open-loop operation.
            obj.AssertConnected();
            obj.device.SetPositionAs(obj.ch1, 0);
            obj.device.SetPositionAs(obj.ch2, 0);
            obj.pos = struct('x', 0, 'y', 0);
        end


        function Jog(obj, axis, direction)
            % One physical jog on the given axis, using whichever raw
            % Jog Size is currently configured on the controller.
            % Confirmed mapping: CH1 Inc=+X, CH1 Dec=-X, CH2 Inc=-Y, CH2 Dec=+Y
            arguments
                obj
                axis (1,1) char {mustBeMember(axis, {'x', 'y'})}
                direction (1,1) double {mustBeMember(direction, [-1, 1])}
            end
            obj.AssertConnected();

            if axis == 'x'
                if direction > 0
                    obj.device.Jog(obj.ch1, obj.jogIncrease, obj.timeout_ms);
                else
                    obj.device.Jog(obj.ch1, obj.jogDecrease, obj.timeout_ms);
                end
                obj.pos.x = obj.pos.x + direction;
            else
                if direction > 0
                    obj.device.Jog(obj.ch2, obj.jogDecrease, obj.timeout_ms);
                else
                    obj.device.Jog(obj.ch2, obj.jogIncrease, obj.timeout_ms);
                end
                obj.pos.y = obj.pos.y + direction;
            end

            pause(obj.settle_s);
        end


        function MoveByDistance(obj, axis, distance_um, stepUm)
            % Moves ~distance_um along axis using repeated Jog() calls,
            % assuming stepUm per jog (the raw Jog Size currently
            % configured on the controller must actually correspond to
            % stepUm -- this function does not verify that). Rounds to
            % the nearest whole jog, so actual distance can be off by up
            % to +/- stepUm/2.
            arguments
                obj
                axis (1,1) char {mustBeMember(axis, {'x', 'y'})}
                distance_um (1,1) double
                stepUm (1,1) double {mustBePositive}
            end
            if distance_um == 0, return; end
            direction = sign(distance_um);
            numJogs = round(abs(distance_um) / stepUm);
            for k = 1:numJogs
                obj.Jog(axis, direction);
            end
        end


        function ReturnToZero(obj)
            % Moves both axes back to the software-tracked (0,0), using
            % whatever Jog Size is currently configured. Does not
            % guarantee a return to the true physical starting point --
            % see the class-level note on open-loop operation.
            try
                while obj.pos.x > 0, obj.Jog('x', -1); end
                while obj.pos.x < 0, obj.Jog('x', +1); end
                while obj.pos.y > 0, obj.Jog('y', -1); end
                while obj.pos.y < 0, obj.Jog('y', +1); end
            catch ME
                warning('KIM101Motor:ReturnFailed', 'Could not return to zero: %s', ME.message);
            end
        end


        function posStruct = GetPosition(obj)
            % Returns the SOFTWARE step counter (obj.pos), NOT a physical
            % measurement -- see the class-level note.
            posStruct = obj.pos;
        end


        function counterCh1Ch2 = GetControllerPosition(obj)
            % BEST-EFFORT: reads the controller's OWN internal step
            % counter for both channels, if the installed Kinesis API
            % exposes it. This can catch a MATLAB<->controller command
            % mismatch (e.g. a dropped Jog command), but it CANNOT catch
            % actual mechanical step loss/slip -- the controller has no
            % more physical feedback than we do.
            obj.AssertConnected();
            try
                c1 = double(obj.device.GetPosition(obj.ch1));
                c2 = double(obj.device.GetPosition(obj.ch2));
                counterCh1Ch2 = [c1, c2];
            catch ME
                error('KIM101Motor:GetPositionUnsupported', ...
                    ['GetPosition is not available on this Kinesis API version (%s). ' ...
                     'Run methods(obj.device) to find the correct call.'], ME.message);
            end
        end


        function SetJogStepSizeUm(obj, stepUm)
            % BEST-EFFORT: sets the jog step size (raw controller steps,
            % converted from stepUm via obj.umPerStep) on both channels
            % via the Kinesis .NET API, instead of requiring manual
            % reconfiguration in the Kinesis GUI.
            %
            % NOTE: exact method/property names for jog parameters can
            % vary between Kinesis .NET API versions. If this throws, run
            %   methods(obj.device)
            % in MATLAB and look for the Get/Set jog-parameter calls in
            % your installed version, then adjust the two GetJogParams/
            % SetJogParams lines below. Callers should wrap this in
            % try/catch and fall back to whatever step size is already
            % configured on the controller if it fails.
            arguments
                obj
                stepUm (1,1) double {mustBePositive}
            end
            obj.AssertConnected();

            rawSteps = round(stepUm / obj.umPerStep);

            jogParamsCh1 = obj.device.GetJogParams(obj.ch1);
            jogParamsCh1.StepSize = rawSteps;
            obj.device.SetJogParams(obj.ch1, jogParamsCh1);

            jogParamsCh2 = obj.device.GetJogParams(obj.ch2);
            jogParamsCh2.StepSize = rawSteps;
            obj.device.SetJogParams(obj.ch2, jogParamsCh2);

            pause(0.1);
        end


        function Close(obj)
            if ~obj.isConnected, return; end
            try, obj.device.StopPolling(); catch, end
            try
                obj.device.Disconnect();
            catch
                try, obj.device.ShutDown(); catch, end
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

        function AssertConnected(obj)
            if ~obj.isConnected
                error('KIM101Motor:NotConnected', 'The KIM101 is not connected.');
            end
        end

    end
end
