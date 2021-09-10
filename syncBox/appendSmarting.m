function allData = appendSmarting(allAcc, window)
% appendSmarting() - append EEG data from a Smarting amp to already
%                   synchronized Faros accelerometers by finding pulses sent by the syncBox 
%                   with a template matching procedureif EEG data only has 1 synchronization
%                   trigger, but Faros have multiple and can be warped.
%           
% v1.0,july 2019
% 
%--------------------------
% Attention: beta version!
%--------------------------
%
% Usage: >> allData = pp_sync(allAcc, window);
%
% Inputs:
%   allAcc       - a cell with EEGLAB structures of EEG & accelerometer data with fields
%       data        - continous, chan x frames,
%       srate       - sampling rate,
%       TTL.chan    - index of channel on which the TTL pulse was recorded,
%       TTL.signal['Smarting'|'LiveAmp'|'accelerometer'] - type of stream
%   window          - timie window in s in which trigger was sent
% 
% Outputs:
%   allData - EEGLAB structure of synchronized data. TTL chan of accelerometers will be deleted,
%            resampled to EEG frequency
%
% naj, April 2019
%
%
% Caveats:
%----------------------------------------------------------------------------------
% Failed template matching:
% - Signal is to noisy (check data before submitting to function)

if nargin < 1
    error('Provide EEG and accelerometer data'); end

% get folder where function is stored
PATHIN = mfilename('fullpath');
PATHIN = PATHIN(1:end-14);
% PATHIN = 'O:\projects\all_gait\Scripts\syncBox\';

%% check input and set sampling frequency for synchronized streams
numAcc = length(allAcc); sratetmp = nan(length(allAcc),1);idxEEG = [];

for s = 1:numAcc
    Acc = allAcc{s};
    % check input structures
    if ~isfield(Acc, 'TTL')
        error(['Structure ', num2str(2), ' does not contain the field "TTL"!']);end
    if ~isfield(Acc.TTL, 'chan')
        error(['Structure ', num2str(2), ' does not contain the field "TTL.chan"!']);end
    if ~isfield(Acc.TTL, 'signal')
        error(['Structure ', num2str(2), ' does not contain the field "TTL.signal"!']);end
    
    % get all sampling rates
    sratetmp(s) = Acc.srate;
    
    % set sampling rate to EEG sampling rate if one of the inputs is EEG
    if strcmp(Acc.TTL. signal, 'Smarting') %if there are multiple EEG datasets the last one is chosen
        idxEEG(end+1) = s;
        srate = Acc.srate;
        idxSmart = s;
    end
end

% otherwise use lowest sampling rate
if ~exist('srate', 'var')
    srate = min(sratetmp);
end


%% find sync pulse of each stream
durACC = nan(1,numAcc);

for s = 1:numAcc % loop through all structures
    
    Acc = allAcc{s};
    
    % flag whether template matching will be skipped (for LiveAmp)
    skip = 0;
    % or manually corrected recodings (trigger and DC saturation overlap)
    if ~isempty(Acc.event) && sum(strcmp({Acc.event.type},'sync1'))>=2
        skip = 1;
        disp([Acc.setname, ' already has sync triggers, skipping this dataset']);
    end
    
    % load template
    if strcmp(Acc.TTL.signal, 'Smarting')
        load([PATHIN, 'tmpEEG'],'tmpAcc');
    elseif strcmp(Acc.TTL.signal, 'LiveAmp')
        skip = 1;
    elseif strcmp(Acc.TTL.signal, 'accelerometer')
        load([PATHIN, 'tmpAcc'],'tmpAcc');
    else
        error('Your TTL.signal is not of type "Smarting" , "LiveAmp" or " accelerometer".');
    end
    
    if ~skip
        
        % check whether sampling rate of template and signal match the determined srate
        if Acc.srate ~= srate
            disp(['Resample ', Acc.TTL.signal,' signal...']);
            setname = Acc.setname;
            evalc('Acc = pop_resample(Acc, srate);');% resample template
            Acc.setname = setname;
        end
        
        if tmpAcc.srate ~= srate
            disp('Resample sync-pulse template...');
            evalc('tmpAcc = pop_resample(tmpAcc, srate);');% resample template
        end
        sigtmpAcc = tmpAcc.data;
        
        
        % get data from channel with TTL pulse & specified time window
        win = window*Acc.srate;
        sigAcc = Acc.data(Acc.TTL.chan,win(1):win(2));
        
        % template matching
        % calculate correlation of template and signal at each frame
        iterations = length(sigAcc)-length(sigtmpAcc);
        rACC = nan(1, iterations);
        disp('Correlating window and template...')
        for i = 1:iterations
            rACCtmp = corrcoef(sigAcc(1,i:i+length(sigtmpAcc)-1),sigtmpAcc);
            rACC(i) = rACCtmp(1,2);
        end
        
        [~, lat] = sort(rACC, 'descend');% get highest correlation
        % look for steepest incline @+- 10 samples
        idx = lat(1)-10;
        [~, iI] = max(diff(sigAcc(idx:idx+20)));% get steepest incline
        Acc.event(end+1).latency = win(1)+idx+iI; %save latency
        Acc.event(end).type = 'sync1';
        Acc.event(end).duration = 1;
        Acc.event(end).correlation = rACC(idx+iI); % ...and their correlation with the template
        
        disp('done')
    end
    
    
    % resample to set srate
    if srate~=Acc.srate
        setname = Acc.setname;
        evalc('Acc = pop_resample(Acc,srate);');
        Acc.setname = setname;
    end


% duration of recording (samples)
evLat = geteventlat(Acc, 'sync1');
durACC(s) = evLat(end)-evLat(1);
Acc.TTL.dur = durACC(s);
if durACC(s) == 0
    idx = find(strcmp({Acc.event.type}, 'sync1'),1); % index of fist sync event
    tmp = Acc.event(idx).latency; % its latency
    durACC(s) = Acc.pnts-tmp;
end

evalc('Acc = eeg_checkset(Acc,"eventconsistency");');
allAcc{s} = Acc;

end


%% find differences in recording duration of EEG and accelerometers
% if exist('idxEEG', 'var')
%     offset = durACC-durACC(idxEEG(1));
%     % if no EEG data is available, offset to first signal is calculated
% else
%     offset = durACC-durACC(1);
% end


%% determine legth of recording
% choose longest one do all data ist there (for warping)
[dur, ~] = min(durACC); %before max


%% prune
for a = 1:numAcc
    % cut signal from first to last (including the pulse itself) TTL pulse and delete TTL pulse channel
    idx = find(strcmp({allAcc{a}.event.type}, 'sync1'),1); % index of fist sync event
    tmp = allAcc{a}.event(idx).latency; % its latency
    evalc('allAcc{a} = pop_select( allAcc{a},"point",[tmp-1 tmp+dur])'); % cut signal
end

%% plot
    figure()
    set(gcf, 'Position', [1 1 1500 800]);
    leg = cell(1, numAcc);
    for a = 1:numAcc
        Acc = allAcc{a};
        data = Acc.data(Acc.TTL.chan,:);
        subplot(211)
        plot(data); hold on
        
        subplot(212)
        plot(data); hold on
        leg{a} = allAcc{a}.setname;
    end
    
    subplot(211)
    xlim([-10, length(sigtmpAcc)+10]);
    xlabel('samples')
    title({'Alignment of streams after synchronization';''; 'Beginning'})
    
    subplot(2, 1, 2)
    xlim([length(allAcc{a}.data)-length(sigtmpAcc)-10, length(allAcc{a}.data+10)]);
    xlabel('samples')
    title('End')
    
    legend(leg,'location', [0.87 0.5 0.12 0.06]);


%% merge
disp('Merging datastreams')
allAcctmp = cell(1,numAcc);
for a = 1:numAcc
    
    % delete TTL chan
    if ~ismember(a,idxEEG)
        evalc('allAcctmp{a} = pop_select( allAcc{a}, "nochannel", allAcc{a}.TTL.chan);');
    else
        allAcctmp{a}= allAcc{a};
    end
    
    % delete not shared event fields
    try allAcctmp{a}.event= rmfield( allAcctmp{a}.event,{'channel','bvtime','bvmknum','code','urevent'}); end
    try allAcctmp{a}.event = rmfield( allAcctmp{a}.event,'correlation');end
    try allAcctmp{a}.event = rmfield( allAcctmp{a}.event,'urevent');end
    
    % and not sheared chalocs info
    try allAcctmp{a}.chanlocs= rmfield( allAcctmp{a}.chanlocs,{'sph_phi_besa','sph_theta_besa'});end
    
    % add channel source (if not present already)
    % not tested yet!
    if ~isfield(allAcctmp{a}.chanlocs, 'source')
        [allAcctmp{a}.chanlocs(:).source] = deal(allAcctmp{a}.setname);
    end
    

    % merge
    if a == 1
        allData = allAcctmp{a};
    else
        % add channel names
        allData.chanlocs = [allData.chanlocs, allAcctmp{a}.chanlocs];
        
        % append data
        allData.data = [allData.data; allAcctmp{a}.data];
        
        % append events
        allData.event = [allData.event, allAcctmp{a}.event];
    end
end

% put in correct order & create event.urevent
evalc('allData = eeg_checkset(allData,"eventconsistency");');
% create index
urevent = num2cell(1:length(allData.event));
% add
[allData.event(:).urevent] = deal(urevent{:});

% adapt channel number...
allData.nbchan = size(allData.data, 1);

% and setname
allData.setname = 'Merged data';
disp('done!');

end