%% F3_jos_sync2_PlotERSP5_generate_brainmeasures2
% find the spatial and frequency ROIs and generate matrices for R analysis

clear all
close all

EXPERIMENT = 'jos_sync';
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'data\ana3_processed\ERSP4\'];
PATHOUT = [MAINPATH, 'data\ana3_processed\data_measures\'];

ETYPE = {'act', 'pas'};
HUMANS = {'Par', 'Exp'};
CONDS_long = {'natural','control', 'sync'};
CNDS =  {'WN','WC','WS'};
CONDS = {'HS_WN','HS_WC', 'HS_WS' };
allGaitEV  = {'TO', 'HS'}; %stepping info w/o conds
allGaitEV2  = {'HS','TO', 'HS'};

RESPONSEtime = [50 2000]; % define time window in which resonse (aka next RHS) can occur [ms]
newLat      = [1 68 100];
newLat2     = [0 680 1000];
load([PATHIN, 'TF5_680both_Par'])

% next steps:
% for each subject:
% Db convert both baseline and tfdata
% average all 3 cond baselines together
% remove averaged baselines from each cond FOR EACH SUBJECT
% plot ERSP

% compute baseline
% baseline:(ch,freqs,i_sub,i_cond)
% data: (ch,freqs,i_sub,times,i_cond)
% db_base (ch, freqs, humans)

% baseline humans (step epochs related to par and exp steps) together
for i_sub = 1:length(TF.subjectNames)
    for ch = 1:length(TF.chanlocs);
        db_base(ch, :, i_sub) =    nanmean([nanmean(TF.baseline(ch, :, i_sub,1, :),5);   nanmean(TF.baseline(ch, :, i_sub,2, :),5);   nanmean(TF.baseline(ch, :, i_sub,3, :),5) ]); % extract data freqs x times for each sub
        
        % baseline correct
        for i_cond = 1:length(CNDS);
            db_tfdata_baseremove(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TF.data(ch, :,i_sub,:,i_cond)),  squeeze(db_base(ch, :, i_sub))'));
        end
    end
end
% stack all conds together
db_tfdatab_combconds = mean(mean([db_tfdata_baseremove(:, :,:,:, :)  ],4),5);
db_tfdatab_combcondsabs = mean(abs(mean(db_tfdata_baseremove(:, :,:,:, :),4))  ,5);

db_tfdatab_combsubs = mean([db_tfdata_baseremove(:, :,:,:, :)  ],5);
db_tfdatab_combsubsabs = abs(mean(db_tfdata_baseremove(:, :,:,:, :)  ,5));

% locs for topo ERSP plot
loc =[ nan nan nan nan   1  33  32 nan nan nan nan ...
    nan nan 5  35     34 nan 63 64 31 nan nan ...
    nan nan 4 36  3  2  29 62 30 nan nan ...
    6  37 nan 38 7  nan 28 60 nan 61 27 ...
    nan 9  40 8  39 24 59 25 58 26 nan ...
    10 41 11 42 12 55 23 56 22 57 21 ...
    44 15 43 14 16 13 18 19 54 20 53 ...
    nan 45 46 47 48 17 49 50 51 52 nan  ];

%% Cluster-averaged ERSPs
% 1. ERSP for all trials
% figure params
lab = {'RHS', 'RTO', 'RHS'};
fnsize = 5;
lat = TF.events.latency;
for i_fig = 1:2;
    % i_human = 1;
    figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
    for i_chan = 1:length(loc)
        if isnan(loc(i_chan));
        else
            if i_fig ==1;
                tmp = squeeze(db_tfdatab_combconds(loc(i_chan), :,:));
                tmp2 = squeeze(mean(db_tfdatab_combconds(:, :,:))); % average of all channel activity
                clim = .3;
            else
                tmp = squeeze(db_tfdatab_combcondsabs(loc(i_chan), :,:));
                tmp2 = squeeze(mean(db_tfdatab_combcondsabs(:, :,:))); % average of all channel activity
                clim = .5;
            end
            subplot(8,11,i_chan)
            contourf(TF.times, TF.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
            plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
            
            caxis([-clim clim]);
            colormap  jet, box off
            L3 =     line([get(gca, 'Xlim')],[8 8],'color','k');
            L4 =     line([get(gca, 'Xlim')],[12 12],'color','k');
            %xticks([]); yticks([])
            yticks(10:10:60);ylabel('Frequency (Hz)');
            title([ TF.chanlocs(loc(i_chan)).labels ])
            set(gca,'FontSize',fnsize)
        end
    end
    
    subplot(8,11,1)
    contourf(TF.times, TF.freqs,tmp2, 50,'linecolor','none'); hold on % plot ERSP
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    
    caxis([-clim clim]);
    colormap  jet, box off
    L3 =     line([get(gca, 'Xlim')],[8 8],'color','k');
    L4 =     line([get(gca, 'Xlim')],[12 12],'color','k');
    %xticks([]); yticks([])
    yticks(10:10:60);ylabel('Frequency (Hz)');
    title([ 'All channels mean' ])
    set(gca,'FontSize',fnsize)
    xticks(lat); xticklabels(lab);xlabel('Time')
    
    caxis([-clim clim]);
    h = colorbar('Position',[0.9 0.44 0.04 0.09], 'FontSize', 8);%('Location', 'east');%('Position',[0.9 0.44 0.04 0.09]);
    ylabel(h, 'Power (a.u.)')
    
end

%% 2. Beta cluster
% average and abs average

% figure params
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
fnsize = 10;
cluschans = [39 24 59];
part = {'mean', 'abs mean'};

for i_fig = 1:2
    figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
    for i_chan = 1:length(cluschans)
        if i_fig == 1;
            tmp = squeeze(db_tfdatab_combconds(cluschans(i_chan), :,:));
            clim = .5;
        else
            tmp = squeeze(db_tfdatab_combcondsabs(cluschans(i_chan), :,:));
            clim = .5;
        end
        subplot(2,3,i_chan)
        contourf(TF.times, TF.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
        plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
        
        caxis([-clim clim]);
        colormap  jet, box off
        L3 =     line([get(gca, 'Xlim')],[16 16],'color','k');
        L4 =     line([get(gca, 'Xlim')],[32 32],'color','k');
        L3 =     line([0 0],[get(gca, 'Ylim')],'color','k');
        yticks(10:10:60);ylabel('Frequency (Hz)');
        title([ TF.chanlocs(cluschans(i_chan)).labels ])
        set(gca,'FontSize',fnsize)
    end
    % average of all channel activity
    if i_fig == 1;
        tmp = squeeze(mean(db_tfdatab_combconds(cluschans, :,:)));
    else
        tmp = squeeze(mean(db_tfdatab_combcondsabs(cluschans, :,:)));
    end
    subplot(2,3,5)
    contourf(TF.times, TF.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    
    caxis([-clim clim]);
    colormap  jet, box off
    L3 =     line([get(gca, 'Xlim')],[16 16],'color','k');
    L4 =     line([get(gca, 'Xlim')],[32 32],'color','k');
    L3 =     line([0 0], [get(gca, 'Ylim')],'color','k');
    yticks(10:10:60);ylabel('Frequency (Hz)');
    title([ 'Cluster ', part{i_fig} ])
    set(gca,'FontSize',fnsize)
    xticks(lat); xticklabels(lab);xlabel('Time')
    
    caxis([-clim clim]);
    h = colorbar('Position',[0.9 0.44 0.04 0.09]);%('Location', 'east');%('Position',[0.9 0.44 0.04 0.09]);
    ylabel(h, 'Power (a.u.)')
    
    
end
%% find the beta range
cluschans = [39 24 59];

figure
tmp = squeeze(mean(db_tfdatab_combcondsabs(cluschans, :,:)));
plot(TF.freqs, mean(tmp'));
title('abs power over frequency')
ylabel('Power (a.u.)')
xlabel('Frequency (Hz)')

% beta seems to have a burst of activity between 14 and 33 Hz.
[a b] = min(abs(TF.freqs-13.5));
[a c] = min(abs(TF.freqs-32.7));

% Full width half max calculation
tmptmp = mean(tmp');
abstrough = (tmptmp(b)+ tmptmp(c))/2;
[thepeak f] = max(tmptmp(b:c));
halfmax = abstrough + ((thepeak - abstrough)/2);
[d e]= min(abs(tmptmp(b:b+f-1)-halfmax));
hm1 = b+e-1; % -1 for all matrix values added to b, because without it, I am counting the first value (b) twice
[d g]= min(abs(tmptmp(b+20+f-1:c)-halfmax));% tweaked to skip the second peak/trough
hm2 = b+20+f-1+g-1;
figure; plot(TF.freqs, mean(tmp'),'b'); hold on;
plot(TF.freqs(hm1:hm2), tmptmp(hm1:hm2), 'r')
plot([TF.freqs(b) TF.freqs(c)], [tmptmp(b) tmptmp(c)], '*b');
% plot(TF.freqs(c), tmptmp(c), '*b');
plot(TF.freqs(b+f-1), abstrough, '*c');
plot(TF.freqs(b+f-1), thepeak, '*m');
plot(TF.freqs(b+f-1), halfmax, '*k');
plot([TF.freqs(hm1) TF.freqs(hm2)], [tmptmp(hm1) tmptmp(hm2)], 'or');
title('abs power over frequency (FWHM in red)')
ylabel('Power (a.u.)')
xlabel('Frequency (Hz)')
legend({'abs power', 'beta range', 'min values', 'mean min value','peak value', 'halfmax value', 'nearest vals to halfmax'}, 'Location', 'SouthEast')

tmp = squeeze(mean(db_tfdatab_combcondsabs(cluschans, :,:)));
figure;plot(TF.times, mean(tmp(hm1:hm2, :)));
title(['abs beta (', num2str(TF.freqs(hm1)),'-', num2str(TF.freqs(hm2)), ' Hz) over time'])
ylabel('Power (a.u.)')
xlabel('Time (ms)')

tmp = squeeze(mean(db_tfdatab_combconds(cluschans, :,:)));
figure; plot(TF.times, mean(tmp(hm1:hm2, :)));
title(['mean beta (', num2str(TF.freqs(hm1)),'-', num2str(TF.freqs(hm2)), ' Hz) over time'])
ylabel('Power (a.u.)')
xlabel('Time (ms)')
%% 2. Alpha/mu mean Right lateral cluster
% average and abs average

% figure params
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
fnsize = 10;
cluschans = [ 25 58 56 22 ];
part = {'mean', 'abs mean'};

for i_fig = 1:2;
    figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
    for i_chan = 1:length(cluschans)
        if i_fig == 1
            tmp = squeeze(db_tfdatab_combconds(cluschans(i_chan), :,:)); % channel ERSP
            tmp2 = squeeze(mean(db_tfdatab_combconds(cluschans, :,:))); %mean over channels
            clim = .5;
        else
            tmp = squeeze(db_tfdatab_combcondsabs(cluschans(i_chan), :,:));
            tmp2 = squeeze(mean(db_tfdatab_combcondsabs(cluschans, :,:))); %mean over channels
            clim = .5;
        end
        
        plotspot = i_chan;
        if plotspot > 2
            plotspot = plotspot + 1;
        end
        % Channel plots
        subplot(2,3,plotspot)
        contourf(TF.times, TF.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
        plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
        
        caxis([-clim clim]);
        colormap  jet, box off
        L3 =     line([get(gca, 'Xlim')],[7.5 7.5],'color','k');
        L4 =     line([get(gca, 'Xlim')],[12.5 12.5],'color','k');
        L3 =     line([0 0],[get(gca, 'Ylim')],'color','k');
        %xticks([]); yticks([])
        yticks(10:10:60);ylabel('Frequency (Hz)');
        title([ TF.chanlocs(cluschans(i_chan)).labels ])
        set(gca,'FontSize',fnsize)
    end
    
    % average of all channel activity
    subplot(2,3,6)
    contourf(TF.times, TF.freqs,tmp2, 50,'linecolor','none'); hold on % plot ERSP
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    
    caxis([-clim clim]);
    colormap  jet, box off
    L3 =     line([get(gca, 'Xlim')],[7.5 7.5],'color','k');
    L4 =     line([get(gca, 'Xlim')],[12.5 12.5],'color','k');
    L3 =     line([0 0],[get(gca, 'Ylim')],'color','k');
    yticks(10:10:60);ylabel('Frequency (Hz)');
    title([ 'Cluster ', part{i_fig} ])
    set(gca,'FontSize',fnsize)
    xticks(lat); xticklabels(lab);xlabel('Time')
    
    caxis([-clim clim]);
    h = colorbar('Position',[0.9 0.44 0.04 0.09]);%('Location', 'east');%('Position',[0.9 0.44 0.04 0.09]);
    ylabel(h, 'Power (a.u.)')
    
end

%% 3. Alpha/mu mean left lateral cluster
% average and abs average

% figure params
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
fnsize = 10;
cluschans = [40 8 11 42];
part = {'mean', 'abs mean'};

for i_fig = 1:2;
    figure('Renderer', 'painters', 'Position', [10 10 1100 900]);
    for i_chan = 1:length(cluschans)
        if i_fig == 1
            tmp = squeeze(db_tfdatab_combconds(cluschans(i_chan), :,:)); % channel ERSP
            tmp2 = squeeze(mean(db_tfdatab_combconds(cluschans, :,:))); %mean over channels
            clim = .5;
        else
            tmp = squeeze(db_tfdatab_combcondsabs(cluschans(i_chan), :,:));
            tmp2 = squeeze(mean(db_tfdatab_combcondsabs(cluschans, :,:))); %mean over channels
            clim = .5;
        end
        
        plotspot = i_chan;
        if plotspot > 2
            plotspot = plotspot + 1;
        end
        subplot(2,3,plotspot)
        contourf(TF.times, TF.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
        plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
        
        caxis([-clim clim]);
        colormap  jet, box off
        L3 =     line([get(gca, 'Xlim')],[7.5 7.5],'color','k');
        L4 =     line([get(gca, 'Xlim')],[12.5 12.5],'color','k');
        L3 =     line([0 0],[get(gca, 'Ylim')],'color','k');
        %xticks([]); yticks([])
        yticks(10:10:60);ylabel('Frequency (Hz)');
        title([ TF.chanlocs(cluschans(i_chan)).labels ])
        set(gca,'FontSize',fnsize)
    end
    
    % average of all channel activity
    subplot(2,3,6)
    contourf(TF.times, TF.freqs,tmp2, 50,'linecolor','none'); hold on % plot ERSP
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    
    caxis([-clim clim]);
    colormap  jet, box off
    L3 =     line([get(gca, 'Xlim')],[7.5 7.5],'color','k');
    L4 =     line([get(gca, 'Xlim')],[12.5 12.9],'color','k');
    L3 =     line([0 0],[get(gca, 'Ylim')],'color','k');
    yticks(10:10:60);ylabel('Frequency (Hz)');
    title([ 'Cluster ', part{i_fig} ])
    set(gca,'FontSize',fnsize)
    xticks(lat); xticklabels(lab);xlabel('Time')
    
    caxis([-clim clim]);
    h = colorbar('Position',[0.9 0.44 0.04 0.09]);%('Location', 'east');%('Position',[0.9 0.44 0.04 0.09]);
    ylabel(h, 'Power (a.u.)')
    
end
%% % find the alpha range
% find the left and right alpha channels

cluschans = [40 8 11 42 25 58 56 22]; % just average them all together
tmp = squeeze(mean(db_tfdatab_combcondsabs(cluschans, :,:)));

figure; plot(TF.freqs, mean(tmp'));
title('abs power over frequency in right parietal cluster')
ylabel('Power (a.u.)')
xlabel('Frequency (Hz)')
ylim([0 .6])

% alpha seems to have a burst of activity between 5 and 11 Hz.
[a b] = min(abs(TF.freqs-3));
[a c] = min(abs(TF.freqs-15));

% Full width half max calculation
tmptmp = mean(tmp');
abstrough = (tmptmp(b)+ tmptmp(c))/2;
[thepeak f] = max(tmptmp(b:c));
halfmax = abstrough + ((thepeak - abstrough)/2);
[d e]= min(abs(tmptmp(b:b+f-1)-halfmax));
hm1 = b+e-1;  %adjusted for symmetry
[d g]= min(abs(tmptmp(b+f-1:c)-halfmax));
hm2 = b+f-1+g-1;
figure; plot(TF.freqs, mean(tmp'),'b'); hold on;
plot(TF.freqs(hm1:hm2), tmptmp(hm1:hm2), 'r')
plot([TF.freqs(b) TF.freqs(c)], [tmptmp(b) tmptmp(c)], '*b');
% plot(TF.freqs(c), tmptmp(c), '*b');
plot(TF.freqs(b+f-1), abstrough, '*c');
plot(TF.freqs(b+f-1), thepeak, '*m');
plot(TF.freqs(b+f-1), halfmax, '*k');
plot([TF.freqs(hm1) TF.freqs(hm2)], [tmptmp(hm1) tmptmp(hm2)], '*r');
title('abs power over frequency (FWHM in red)')
ylabel('Power (a.u.)')
xlabel('Frequency (Hz)')
legend({'abs power', 'beta range', 'min values', 'mean min value','peak value', 'halfmax value', 'nearest vals to halfmax'})

tmp = squeeze(mean(db_tfdatab_combcondsabs(cluschans, :,:)));
figure;plot(TF.times, mean(tmp(hm1:hm2, :)));
ylim([0 .55])
hold on;   plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
title(['abs alpha ', num2str(TF.freqs(hm1)),'-', num2str(TF.freqs(hm2)), ' Hz) over time'])
ylabel('Power (a.u.)')
xlabel('Time (ms)')

tmp = squeeze(mean(db_tfdatab_combconds(cluschans, :,:)));
figure; plot(TF.times, mean(tmp(hm1:hm2, :)));
ylim([-.25 .25])
hold on;   plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
title(['mean alpha (', num2str(TF.freqs(hm1)),'-', num2str(TF.freqs(hm2)), ' Hz) over time'])
ylabel('Power (a.u.)')
xlabel('Time (ms)')

%% Topos ROI averaged over time
tfranges = [16 32;  7.5 12.5];

figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
for i_col = 1:2;
    for i_band = 1:length(tfranges)
             
        %freq and time range
        range = find(TF.freqs>tfranges(i_band,1),1)-1:find(TF.freqs>tfranges(i_band,2),1)-2;
        timerange = find(TF.times>0,1)-1:find(TF.times>992,1)-2;
        if i_col == 1;
            tmp = mean(mean(db_tfdatab_combconds(:, range, timerange),2),3);
            clim = .1;
        else
            tmp = mean(mean(db_tfdatab_combcondsabs(:, range, timerange),2),3);
            clim = .3;   
        end
        [sth, sorttmp] = sort(abs(tmp));
        bigvals = sorttmp(59:64);
        
        subplot(4,2, 2*(i_band)+i_col) % starts at second row because supertitle warps the top row for some reason
        topoplot(tmp, TF.chanlocs, 'maplimits', [-clim clim],'emarker2',{[bigvals],'s','m'} );
        title([num2str(tfranges(i_band,1)), '-' num2str(tfranges(i_band,2)), 'Hz, ', 'time-averaged ']);
    end
    caxis([-clim clim]);
    h1 = colorbar('Position',[0.5+0.4*(i_col-1) 0.44 0.04 0.09]); % position: [left, bottom, width, height]
    ylabel(h1, 'Power (a.u.)')
end

%% topoplot video
tfranges = [16 32; 7.5 12.5];
lat = [0  .680 1.000];

part = {'mean'};

for i_band = 1:length(tfranges);
    f = figure('Renderer', 'painters', 'Position', [10 10 800 600]);
    %freq and time range
    range = find(TF.freqs>tfranges(i_band,1),1)-1:find(TF.freqs>tfranges(i_band,2),1)-2;
    timerange = find(TF.times>-200)-1;
    tmp = squeeze(mean(db_tfdatab_combconds(:, range, timerange),2));
    clim = .3;
    
    [mymov, colmap] =  eegmovie(tmp, 100, TF.chanlocs,'title', [num2str(tfranges(i_band,1)), '-', num2str(tfranges(i_band,2)), ' Hz'], 'minmax',[-clim, clim],'startsec',-.2,'framenum', 'off', 'time', 'on','vert', lat );
    % seemovie(mymov)
    
    % save movie
    vidObj = VideoWriter(['step_', num2str(tfranges(i_band,1)), '-', num2str(tfranges(i_band,2)), '_Hz.mp4'], 'MPEG-4');
    vidObj.FrameRate = 15;
    open(vidObj);
    writeVideo(vidObj, mymov);
    close(vidObj);
end

%% Histograms: Beta
% data: (ch,freqs,i_sub,times,i_cond)
i_band = [16 32];
cluschans = [39 24 59];

% parameters
frange = find(TF.freqs>i_band(1))-1:find(TF.freqs>i_band(2))-2;
timerange = find(TF.times>0)-1:find(TF.times==1000); %-2 not used cuz 1000 is at the edge of the tf plot (not the epoch)

% graph them averaged
tmp = mean(mean(mean(db_tfdatab_combsubs(cluschans, frange, timerange, :),3),2),1);
figure;
subplot(2,1,1)
histogram(tmp, 15)
title('distribution of subject values (averaged over conditions)')
ylabel('number of subs (out of 18)')
xlabel('beta power averaged over the walking cycle')

% graph them separated
subplot(2,1,2)
tmp1 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans, frange, timerange, :,1),3),2),1));
tmp2 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans, frange, timerange, :,2),3),2),1));
tmp3 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans, frange, timerange, :,3),3),2),1));
histogram(tmp1, 15); hold on;
histogram(tmp2, 15)
histogram(tmp3, 15)
title('distribution of subject values (within conditions)');
ylabel('number of subs (out of 18)');
xlabel('beta power averaged over the walking cycle');
legend({'natural', 'control', 'sync'});

% make the table and save
CONDS_id = {'id', 'a_natural', 'b_control', 'c_sync'};
tab = table( TF.subjectNames, tmp1, tmp2, tmp3, 'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'betacentr.csv']);
betacentr = [tmp1 tmp2 tmp3];

%beta descriptive stats
[h p ci test] = ttest(betacentr(:,1), betacentr(:,2));
Mdiff = mean(betacentr(:,1))-mean(betacentr(:,2));
disp('nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%beta descriptive stats
[h p ci test] = ttest(betacentr(:,3), betacentr(:,1));
Mdiff = mean(betacentr(:,3))-mean(betacentr(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])
%% Histograms: Alpha lateral
% data: (ch,freqs,i_sub,times,i_cond)
i_band = [7.5 12.5];
cluschans1 = [ 25 58 56 22 ];% right
cluschans2 = [40 8 11 42];% left
frange = find(TF.freqs>i_band(1))-1:find(TF.freqs>i_band(2))-2;
timerange = find(TF.times>0)-1:find(TF.times==1000); %-2 not used cuz 1000 is at the edge of the tf plot (not the epoch)

tmp = mean(mean(mean(db_tfdatab_combsubs(cluschans1, frange, timerange, :),3),2),1);
figure;
subplot(2,1,1)
histogram(tmp, 15)
title('RIGHT SIDE distribution of subject values (averaged over conditions)')
ylabel('number of subs (out of 18)')
xlabel('alpha power averaged over the walking cycle')
xlim([-0.15 0.15])

% right side
subplot(2,1,2)
tmp1 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans1, frange, timerange, :,1),3),2),1));
tmp2 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans1, frange, timerange, :,2),3),2),1));
tmp3 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans1, frange, timerange, :,3),3),2),1));
histogram(tmp1, 15); hold on;
histogram(tmp2, 15)
histogram(tmp3, 15)
title('RIGHT distribution of subject values (within conditions)');
ylabel('number of subs (out of 18)');
xlabel('alpha power averaged over the walking cycle');
legend({'natural', 'control', 'sync'});
xlim([-1 2])

%left side
tmp = mean(mean(mean(db_tfdatab_combsubs(cluschans2, frange, timerange, :),3),2),1);
figure;
subplot(2,1,1)
histogram(tmp, 15)
title('LEFT SIDE distribution of subject values (averaged over conditions)')
ylabel('number of subs (out of 18)')
xlabel('alpha power averaged over the walking cycle')
xlim([-0.15 0.15])

% left side
subplot(2,1,2)
tmp4 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans2, frange, timerange, :,1),3),2),1));
tmp5 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans2, frange, timerange, :,2),3),2),1));
tmp6 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans2, frange, timerange, :,3),3),2),1));
histogram(tmp4, 15); hold on;
histogram(tmp5, 15);
histogram(tmp6, 15);
xlim([-1 2]);
title('LEFT distribution of subject values (within conditions)');
ylabel('number of subs (out of 18)');
xlabel('alpha power averaged over the walking cycle');
legend({'natural', 'control', 'sync'});
CONDS_lr = {'id','a_natural_right', 'b_control_right', 'c_sync_right','d_natural_left', 'e_control_left', 'f_sync_left'};

% % table 1 organized for wilcoxon
tab1 = table( TF.subjectNames, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, 'VariableNames', CONDS_lr);
writetable(tab1,[PATHOUT, 'alphalat1.csv']);

% table 2 organized for R
CONDS_id = {'id', 'hemi', 'a_natural', 'b_control', 'c_sync'};
hemi_tm1 = [tmp4;tmp1]; %left over right
hemi_tm2 = [tmp5;tmp2];
hemi_tm3 = [tmp6;tmp3];
hemi_labs = [repmat({'left'},18,1);repmat({'right'},18,1)];
id_hemi = [repmat(TF.subjectNames, 2,1)];
tab = table( id_hemi,hemi_labs,hemi_tm1, hemi_tm2, hemi_tm3, 'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'alphalat2.csv']);

alpharight = [tmp1 tmp2 tmp3];
alphaleft = [tmp4 tmp5 tmp6];

%alpha right
[h p ci test] = ttest(alpharight(:,1), alpharight(:,2));
Mdiff = mean(alpharight(:,1))-mean(alpharight(:,2));
disp('right: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%alpha right
[h p ci test] = ttest(alpharight(:,3), alpharight(:,1));
Mdiff = mean(alpharight(:,3))-mean(alpharight(:,1));
disp('right: sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%alpha left
[h p ci test] = ttest(alphaleft(:,1), alphaleft(:,2));
Mdiff = mean(alphaleft(:,1))-mean(alphaleft(:,2));
disp('left: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%alpha left
[h p ci test] = ttest(alphaleft(:,3), alphaleft(:,1));
Mdiff = mean(alphaleft(:,3))-mean(alphaleft(:,1));
disp('left: sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])