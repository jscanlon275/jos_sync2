%% calculate footprint with Joannas prepared data
% try to calculate gait artifact footprint w/ Joannas data (jos_sync)
% check data first (ERSPs and ERPs)
% rely on assumptions for left foot gait events (as no gait events from 2nd foot
% available)
clc; clear; close all

% WORKSTATION     = 'LAPTOP'; %'PC' 'LAPTOP'; % machine
CALC_FEATURES   = 1; % check whether correct channel location were indext for feature computation
CHECK_FEATURES  = 1; % boxplots of features
STATS_FEATURES  = 1; % check whether features were diffeerent before and after processing
PLOT_FEATURES   = 1; % make radar plot
CHECK_DATA = 0;

proc = {'filtered', 'clean'}; %processing stage (filteres, cleaned)

%% prepare workstation
        addpath \\C:\Users\ebmi2273\jos_sync2\eeglab2020_0 % add eeglab path, if not saved
        PATH = 'C:\Users\ebmi2273\jos_sync2\'; % project directory

addpath([PATH, 'code\gaitEEGfootprint-master\'])
% link to data----------------------------------------------------
PATHIN = [PATH, 'data\ana3_processed\']; %all data
PATHOUT = [PATH, 'data\footprint\'];
PATHOUTacc = [PATHOUT, 'acc', filesep];
if ~exist(PATHOUT, 'dir')
    mkdir(PATHOUT); end

cd(PATH)
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% footprint parameters----------------------------------------------
% epoching
% [ADAPT] event names and order
params.EV           = {'RightHS'};          % [ADAPT] event to epoch around
params.allGaitEV    = {'RightTO', 'RightHS', 'RightHS'}; % correct order of all gait events that will be used for warping
params.newLatms     = [0 180 500 680 1000]; % [ADAPT] new latencies (to warp to), in ms

% channel indices (dim 1 of the time-frequency decomposed data)
% [ADAPT] channel selection basel on you layout (here: custom 64ch layout)
params.B.lateralChanIdx = [1,4,5,6,9,10,15,16,17,18,20,21,26,27,30,31,32,33,35,37,41,44,45,46,48,49,51,52,53,57,61,64]; % index of channels labelled as lateral
params.E.neckChanR      = [49 52 18 51];  % index of channels located over the right side of the neck
params.E.neckChanL      = [45 48 46 16];  % index of channels located over the left side of the neck

%% calculate feature
% ERP = mean(EEGwalkEP.data,3); load this from my data?
if CALC_FEATURES
    
% test other way of calculating:
    for p = 1:2
        if p == 1 %rawdata
            stand = load([PATHIN, 'ERSP4_reref_preAC/','TF5_680both__stbase_preAC2.mat'], 'TF');
            walk = load([PATHIN, 'ERSP4_reref_preAC/', 'TF5_680both_Par_preAC.mat'], 'TF');
        elseif p == 2 % after artifact attenuation
            stand = load([PATHIN, 'ERSP4/', 'TF5_680both__stbase.mat'], 'TF');
            walk = load([PATHIN, 'ERSP4/', 'TF5_680both_Par.mat'], 'TF');
        end
        
        footprint = table('Size',[length(stand.TF.subjectNames), 5],...
            'VariableTypes',{'double','double','double','double','double'},...
            'VariableNames',{'B','C','D','E','F'});
        
        if ~exist([PATH, 'script\gaitEEGfootprint-master\code\footprintParams.mat'], 'file')% only run once, some for both datasets
            % sample indices of the time-frequency decomposed data
            % calculate relevant points from joannas data!
            params.newLatpnts   = find(ismember(walk.TF.times, params.newLatms));     % new latencies in pnts
            params.E.pntsRHS        = params.newLatpnts(1):params.newLatpnts(2);          % double support following right-heel strike
            params.E.pntsLHS        = params.newLatpnts(3):params.newLatpnts(4);         % double support following left-heel strike
            params.D.pntsDouble     = [params.E.pntsRHS, params.E.pntsLHS];
            save([PATH, 'script\functions\gaitEEGfootprint-master\code\footprintParams.mat'], 'params')
        else
            load([PATH, 'script\functions\gaitEEGfootprint-master\code\footprintParams.mat']);
        end
        
        for s = 1:length(stand.TF.subjectNames)
            from = find(walk.TF.times == 0); % RHS1
            to = find(walk.TF.times == 1000); % RHS2
            TFdata = squeeze(mean(walk.TF.data(:,:,s,from:to,:),5));
            % chans x freqs x subjects x pnts x conds (average over conds)
            % --> chans # freqs x pnts (HS to HS)
            TFbaseline = stand.TF.baseline(1:64,:,s); % chans x freqs x subjects
            TFstandBL = TFdata./repmat(TFbaseline,1,1,size(TFdata,3));
                  
            % B) power ratio lateral/medial channels -------------------------
           powLat = sum(TFstandBL(params.B.lateralChanIdx, :,:),'all');
           powMed = sum(TFstandBL, 'all')-powLat;
           feature(1) = powLat/powMed;
            
            % C) correlation across frequencies --------------------------
            feature(2) = Rfreq(TFdata, TFbaseline);
            
            % D) power double support/single supp gait cycle power -------------
            powD = sum(TFstandBL(:,:,params.D.pntsDouble), 'all');
            powS = sum(TFstandBL, 'all')-powD;
            feature(3) = powD/powS;
            
            % E) power at neck electrodes contralateral to HS/ipsi --------------
            powContra = sum(TFstandBL(params.E.neckChanL,:,params.E.pntsRHS),'all')+...
                sum(TFstandBL(params.E.neckChanR,:,params.E.pntsLHS),'all'); 
            powIpsi = sum(TFstandBL(params.E.neckChanL,:,params.E.pntsLHS),'all')+...
                sum(TFstandBL(params.E.neckChanR,:,params.E.pntsRHS),'all'); 
            feature(4) = powContra/powIpsi;
            
            % F) 1-S/W power ratio --------------------------------------
            feature(5) = mean(TFstandBL, 'all');
            
            footprint{s,:} = feature;
        end
        writetable(footprint, [PATHOUT filesep 'new_features_' proc{p}]);
    end
end
 %% % check feature values
if CHECK_FEATURES 
    
    PATHIN = fullfile(PATH, 'data','footprint');
    PATHOUT = fullfile(PATH, 'results','footprint');
    if ~exist(PATHOUT, 'dir')
        mkdir(PATHOUT); end
    
    % set params
    figure; set(gcf, 'units', 'centimeters', 'Position', [0 0 20 10]);
    
    
    for p = 1:length(proc) %processing stage (filteres, cleaned)
        footprint = readtable([PATHIN filesep 'new_features_' proc{p} '.txt']);
        allFootprint.(proc{p}) = footprint;
        
        % plot boxplots of features
        subplot(1,2,p)
        plot([-0.5 7],[0 0], 'k-'); hold on;
        plot(footprint.Variables', 'c*');
        boxplot(footprint.Variables)
        xticklabels({'B','C','D','E','F'});
        xlabel('Feature'); ylabel('value (a.u.)')
        title({proc{p}, ['N = ', num2str(height(footprint))]});
        ylim([-0 2]); xlim([.5 5.5]);
        grid on
    end
    print(fullfile(PATH, 'results','footprint', 'boxplots_all_features'), '-dpng', '-r300');
    close
    
    save(fullfile(PATH, 'results','footprint', 'all_features_allConds'), 'allFootprint');
end

%% statistical evaluation
if STATS_FEATURES  % check whether features were diffeerent before and after processing
    PATHIN = fullfile(PATH, 'results','footprint');
    load(fullfile(PATHIN, 'all_features_allConds'), 'allFootprint');
    
    distFootprint = nan(18,1);
for s = 1:length(distFootprint)
    distFootprint(s) = norm(allFootprint.clean{s,:}-allFootprint.filtered{s,:});
end

footprint = table(distFootprint, 'VariableNames', {'dist'});
writetable(footprint, [PATHIN filesep 'distances']);

% plot
addpath(fullfile(PATH, 'code','RainCloudPlots-master','tutorial_matlab'));
grey = [.6 .6 .6];
axes_font_size = 10;

figure; set(gcf, 'units', 'centimeters', 'Position',[0 0 10 5])
h = raincloud_plot(distFootprint, 'box_on',1,'box_dodge', 1, 'box_dodge_amount',...
.3, 'dot_dodge_amount', .3, 'color', grey, 'cloud_edge_col', grey,'line_width', 1);
ylabel('density (a.u.)'), xlabel('distance (a.u.)')
xlim([0 3])
box off
set(gca, 'FontSize', axes_font_size)
    
end

%% illustrate results
if PLOT_FEATURES   % make radar plot
    %  get radar plot from https://de.mathworks.com/matlabcentral/fileexchange/59561-spider_plot
    % radar plot parameters
    axes_interval = 6;% Axes properties
    axes_precision = 2;
    axes_display = 'one';
    marker_type = 'none';
    axes_font_size = 10;
    label_font_size = 12;
    axes_labels = {'B','C','D','E','F'}; % Axes labels
    fill_option = 'on';
    fill_transparency = 1;
    orange = [.9 .6 0];
    blue = [0 .45 .7];
    lightorange = [251 240 217]/255;
    lightblue = [217 234 244]/255;
    edgeColors = [orange;blue];
    colors = [lightorange;lightblue];
    line_width = 1;
    line_style ={'--','-'};
 
    PATHIN = fullfile(PATH, 'results','footprint');
    load(fullfile(PATHIN, 'all_features_allConds'), 'allFootprint');
    
    %% average data
    P = [mean(allFootprint.filtered.Variables);mean(allFootprint.clean.Variables)];% extract data (obs x vars)
    axes_limits = repmat([0; 3],1,size(P,2));
    
    % Spider plot
    spider_plot_nj(P,...
        'AxesLabels', axes_labels,...
        'AxesInterval', axes_interval,...
        'AxesPrecision', axes_precision,...
        'AxesDisplay', axes_display,...
        'AxesLimits', axes_limits,...
        'FillOption', fill_option,...
        'FillTransparency', fill_transparency,...%        'Color', colors,...
        'LineWidth', line_width,...
        'LineStyle', line_style,...
        'Marker', marker_type,...
        'Color', colors,...
        'edgecolor', edgeColors,...
        'AxesFontSize', axes_font_size,...
        'LabelFontSize', label_font_size);
    legend1 = legend({'PreAC', 'Post-AC'}, 'Fontsize', label_font_size);
    set(legend1, 'Position',[0.45 0.059 0.18 0.06]);
    legend('boxoff')
    set(gcf, 'units', 'centimeters', 'Position',[0 0 10 15])
%     sgtitle('Mean gait artifact footprint')
    savefig([PATHOUT filesep 'artifact_reduction_footprint_radar'])
    %close
end

%% take a look at the data
if CHECK_DATA

% load info
load(fullfile(PATH, 'code', 'chanSubplotID.mat'))
load([PATH, 'code\footprintParams.mat'])

% set params
lab = params.allGaitEV;
lat = params.newLatms;
climAll = [4, 1.5];
fnsz = 12;

PATHIN = fullfile(PATH, 'data','updated_data','ERSP4');
PATHOUT = fullfile(PATH, 'results','footprint');
if ~exist(PATHOUT, 'dir') % crate output dirctory (if neccessary)
    mkdir(PATHOUT); end

for p = 1:2
    if p == 1 %rawdata
        stand = load([PATHIN,filesep, 'TF5_680both__stbase_preAC2.mat'], 'TF');
        walk = load([PATHIN,filesep, 'TF5_680both_Par_preAC.mat'], 'TF');
    elseif p == 2 % after artifact attenuation
        stand = load([PATHIN,filesep, 'TF5_680both__stbase.mat'], 'TF');
        walk = load([PATHIN, filesep,'TF5_680both_Par.mat'], 'TF');
    end
    from = find(walk.TF.times == 0); % RHS1
    to = find(walk.TF.times == 1000); % RHS2
    TFdata = squeeze(mean(walk.TF.data(:,:,:,from:to,:),5));
    % chans x freqs x subjects x pnts x conds (average over conds)
    % --> chans # freqs x pnts (HS to HS)
    TFbaseline = stand.TF.baseline(1:64,:,:); % chans x freqs x subjects
    
    % baseline correct data (dB change to standing baseline)
%     TFdB = 10*log10(TFdata./repmat(TFbaseline,1,1,1,size(TFdata,4)));
TFdB = (TFdata./repmat(TFbaseline,1,1,1,size(TFdata,4)));
    
    freqs = walk.TF.freqs;
    times = walk.TF.times(from:to);
    
    %%% plot average across subjects %%%
    fig = figure; set(gcf, 'units', 'centimeters','Position', [0 0 40 25])
    set(0, 'DefaultAxesFontSize', 6)
    for ch = 1:length(walk.TF.chanlocs) % plot each chan seperatly
        idxChan = find(strcmp(plotID.ID, walk.TF.chanlocs(ch).labels));
        subplot(10,11,plotID.num(idxChan))
        contourf(times, freqs,squeeze(mean(TFdB(ch,:,:,:),3)), 50,'linecolor','none'); hold on % plot ERSP
        plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
        
        caxis([-climAll(p) climAll(p)]); %        colorbar
        colormap  jet, box off
        
        title(walk.TF.chanlocs(ch).labels)
        xticks([]); yticks([]); yticks(10:20:60);
        set(gca,'FontSize',6)
    end
    subplot(10,11,110)
    contourf(times, freqs,squeeze(mean(mean(TFdB,1),3)), 50,'linecolor','none'); hold on % plot ERSP
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    
    caxis([-climAll(p) climAll(p)]);
    colormap  jet, box off
    
    title('Channel mean')
    xticks(lat);xticklabels({'RHS', 'RTO', 'RHS'});
    yticks([]); yticks(10:20:60);
    set(gca,'FontSize',6)
    
    h = colorbar('Position',[0.92 0.54 0.01 0.05]);
    ylabel(h, 'Power change to standing BL (dB)')
    
    han=axes(fig,'visible','off');
    han.XLabel.Visible='on';xlabel(han,'Time', 'FontSize', fnsz);
    han.YLabel.Visible='on';ylabel(han,'Frequency (Hz)', 'FontSize', fnsz);
    
    sgtitle(['gait ERSPs of all subjects, ', proc{p},' data'],...
        'FontSize', fnsz,'interp', 'none')
    
    print([PATHOUT filesep 'jos_sync_all_ERSPrel_ eachChan_', proc{p}],'-dtiff','-r300' )
%     savefig([PATHOUT filesep 'jos_sync_all_ERSP_ eachChan_', proc{p}])
    close
end

end
