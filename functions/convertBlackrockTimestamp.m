function [timestampSamples, info] = convertBlackrockTimestamp(rawTimestamp, TimeRes, SampleRes)

% AUTHOR: Niloufar Shahdoust; niloufar.shahdoust@utah.edu

% this code converts Blackrock timestamps to old sample indices.

% this code has been written based on Blackrock explanation here:
% https://support.blackrockneurotech.com/portal/en/kb/articles/ptp-alignment:

% which is:
% Timestamping has fundamentally changed from operating on a nominal 30 kHz
% clock to Precision Time Protocol (PTP) on Gemini hardware. This is a 
% departure from previous communication between the digital hubs and the I/O 
% signal processor, where a single clock driven by the Neural Signal 
% Processor (NSP) would sample from the data streaming from the digital 
% hubs via fiber optic, which resulted in data that was properly aligned by default.
% Now, the Gemini hubs and I/O NSP each have their own clocks. 
% This has the benefit of being able to run independently from the I/O module 
% as needed but creates a necessity for ensuring that the data in recorded
% files is aligned with time.


% input:
%   rawTimestamp        - timestamp vector from NEV/NSx/eventTimes
%   TimeRes             - NEV.MetaTags.TimeRes
%   SampleRes           - NEV.MetaTags.SampleRes, my case: 30K

%
% output:
%   timestampSamples    - timestamps converted to sample indices
%   info                - structure with conversion info
%
% example in my code:
%   [spikeTS, info] = convertBlackrockTimestamp( ...
%       NEV.Data.Spikes.TimeStamp, ...
%       NEV.MetaTags.TimeRes, ...
%       NEV.MetaTags.SampleRes);



    TimeRes   = double(TimeRes);
    SampleRes = double(SampleRes);

    % find PTP timestamps
    isPTP = TimeRes ~= SampleRes || max(double(rawTimestamp(:))) > 1e9;

    info = struct();
    info.TimeRes = TimeRes;
    info.SampleRes = SampleRes;
    info.isPTP = isPTP;

    if isPTP

        % Gemini PTP timestamps.
        % they are laaaarrrge numbers so turn to uint64
        rawTimestamp_u64 = uint64(rawTimestamp);

        % the first/minimum timestamp as time zero. because 
        % the new PTP timestamp is like an absolute clock time,
        %  not a sample number from the start of recording.
        referenceTimestamp = min(rawTimestamp_u64(:));

        info.referenceTimestamp = referenceTimestamp;

        % convert PTP time to 30 kHz sample index
        % subtracts the reference timestamp, so the huge PTP timestamp becomes elapsed time from the start
        % multiplies by the SampleRes
        % divides by the TimeRes to convert time units into sample indices
        timestampSamples = round(double(rawTimestamp_u64 - referenceTimestamp) .* SampleRes ./ TimeRes);

    else

        % for old Blackrock files:
        % timestamps are already sample indices.
        timestampSamples = double(rawTimestamp);
        info.referenceTimestamp = 0;

    end
end