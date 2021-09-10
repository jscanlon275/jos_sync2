function [EEG] = naj_steptrig(allDATA, varargin)
% pp_steptrig() - inserts left and right step triggers into EEG streams
% 
% Not optional: Use only after pp_sync function (syncing the accelerometers and EEG)
% 
% Strong recomendation before use: epoch the data so that it approximatelly 
% includes only walking periods. The step detection is robust and steps 
% may be detected outside of the designated task period. If there are large
% artifacts in accelerometer data some steps may not be detected. In any case use 
% plotting option to check visually for misdetected steps or steps outside of
% the walking period.
% 
% Usage: >> [allDATA, AMP, ACC] = pp_sync(ACC, , , plot);
%        >> [allDATA, AMP, ACC] = pp_sync(ACC);
%
% Inputs:
%   ACC      - cell with EEG and accelerometer strructures with fields. In
%              accelerometer data it is assumed that second row is x axis 
%              (back to front).
%       
% Optional inputs:
%     'defaulttemplate'  - Use a predifined step template (value 1) or define a 
%                 template from the current dataset (value 0). Default is 1. 
%     'threshold' - Amplitude threshold of scaled correlation signal. Default is 1.2.
%     'minStep'   - Minimal distance between two crosscorellation peak in seconds.
%                 Default is 0.6.
%     'plottrig'  - Plot the accelerometer signal with the resulting 
%                 triggers (value 1). Default 1.
% 
% Outputs:
%     EEG - EEG structure with step triggers
%     allACC - Accelerometer and EEG structures in the same cell with triggers included.
%     
%
% Bojana Mirkovic, Dec 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% assign inputs

%check which dataset is which and rename
% will throw error if Acc is resampled!!! find better way!!!, supress
% output with evalc()
AccL=pop_select(allDATA,'channel',find(strcmpi({allDATA.chanlocs.source}, 'Accelerometer left')));
AccR=pop_select(allDATA,'channel',find(strcmpi({allDATA.chanlocs.source}, 'Accelerometer right')));
EEG = allDATA;
% EEG=pop_select(allDATA,'channel',find(strcmpi({allDATA.chanlocs.source}, 'EEGresampled')));

UseTemp=[];
PlotTrig=[];
ThresholdConst=[];
Distance=[];
CheckTemp = [];

if nargin < 1
    error('Provide EEG and accelerometer data'); end
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'defaulttemplate') 
        UseTemp= varargin{i+1};
    elseif strcmpi(varargin{i},'plottrig')
        PlotTrig=varargin{i+1}; 
    elseif strcmpi(varargin{i},'threshold')
        ThresholdConst=varargin{i+1}; 
    elseif strcmpi(varargin{i},'minStep')
        Distance=varargin{i+1}*AccL.srate; 
    elseif strcmpi(varargin{i},'CheckTemp')
        Distance=varargin{i+1}; 
    end
end

if isempty(UseTemp)
    UseTemp = 1; end
if isempty(CheckTemp)
    CheckTemp = 0; end
if isempty(PlotTrig)
    PlotTrig = 1; end
if isempty(ThresholdConst)
    ThresholdConst = 1.2; end
if isempty(Distance)
    Distance = 0.6*AccL.srate; end


Path=mfilename('fullpath'); Path=Path(1:end-12);

%% Find or load template

if UseTemp
    load([Path 'DefStepTmp.mat']);

    if (Step.srate~=AccL.srate) || (Step.srate~=AccR.srate)
        warning('Sampling frequency of the loaded template does not match sampling frequency of the accelerometer signal. Resample template')
        Step = pop_resample(Step, AccR.srate);
    end
    
    if CheckTemp
    figure(213)
    subplot(1,4,1:3)
        plot(AccL.times,AccL.data(1,:)); hold on;
        plot(AccL.times,AccR.data(1,:));
        legend('Left','Right')
        title('Zoom in and compare AccX and template; to continue press enter in command window.');
    subplot(1,4,4)
        plot(Step.times, Step.data); hold on;
        title('Chosen template');    
    disp('Look at the figure. Zoom in and compare x axis of accelerometer data and chosen template.')
    disp('To continue to accept/decline template dialog press enter in command window.')
    disp('In case you do not wish to proceed with current template you will be given an option to select a template from current accelerometer data.')
    pause
    %ask if the template fits
    answer = questdlg('Proceed with this template?', 'Template', ...
                       'Yes, proceed.','No, start template selection.', 'Yes, proceed.');
    close(213)

    %save template or restart selection
    if strcmp(answer, 'No, start template selection.')
        UseTemp=0;
    end
    
    end

end
    
if ~UseTemp
    startTemplate=1;
    while startTemplate
        % show plot of Left and right acc data
        figure(212)
            plot(AccL.data(1,:)); hold on;
            plot(AccR.data(1,:))
            legend('Left','Right')
            title('Zoom in to ~3 steps and press enter in command window.')
        pause; 
        title('Choose starting point of the step on left accelerometer')
        [x1,~]=ginput(1); x1=round(x1);
        title('Choose end point of the same step on left accelerometer')
        [x2,~]=ginput(1); x2=round(x2);
        close 212
        
        Step.data=AccL.data(1,x1:x2);
        Step.times=AccL.times(x1:x2)-AccL.times(x1);
        Step.srate=AccL.srate;

        %show chosen template
        figure(212)
            plot(Step.Time/1000,Step.Tmp); 
            xlabel('Time (s)');
            ylabel('Amplitude');
            title('Template step');

        %ask if the template fits
        answer = questdlg('Is this your template?', 'Template', ...
                'Yes, save and proceed.', 'Yes, proceed without saving.', 'No, restart template selection.', 'Yes, proceed without saving.');
        close(212)

        %save template or restart selection
        if strcmp(answer, 'Yes, save and proceed.')
            save([Path 'CurrentStepTemp.mat'],'Step');
            startTemplate=0;
        elseif strcmp(answer, 'Yes, proceed without saving.')
            startTemplate=0;
        elseif strcmp(answer, 'No, restart template selection.')
            startTemplate=1;
        end
    end
            
        
end   

%% find steps in accelerometer data
%length of the template, necessary for crosscorr later
Wdw = length(Step.data); 
disp('Finding steps, please wait. This may take a couple of minutes.')

% Crosscorelation and step detection
for i=1:length(AccL.data(1,:))-Wdw
    AccL.ccorr(i) = xcorr(AccL.data(1,i:i+Wdw), Step.data,0);
    AccR.ccorr(i) = xcorr(AccR.data(1,i:i+Wdw), Step.data,0);
end
AccL.ccorr=bom_scaledata(AccL.ccorr,-100,2000); %function attached, only scales data to a given range
AccR.ccorr=bom_scaledata(AccR.ccorr,-100,2000);
Threshold_L=ThresholdConst*rms(AccL.ccorr);
Threshold_R=ThresholdConst*rms(AccR.ccorr);
[AccL.peakAmp, AccL.peakSmp] = findpeaks(AccL.ccorr,'MinPeakHeight', Threshold_L,'MinPeakDistance',Distance);
[AccR.peakAmp, AccR.peakSmp] = findpeaks(AccR.ccorr,'MinPeakHeight', Threshold_R,'MinPeakDistance',Distance);

%show crosscorr pattern and raw data on the same plot
if PlotTrig
figure(212)
subplot(2,1,1)
    plot(AccL.ccorr); hold on;
    plot(AccL.data(1,:));
    plot(AccL.peakSmp,AccL.peakAmp,'r*')
    title('Left foot acc');
    xlabel('Time (samples)');
subplot(2,1,2)
    plot(AccR.ccorr); hold on;
    plot(AccR.data(1,:));
    plot(AccR.peakSmp,AccR.peakAmp,'r*')
    title('Right foot acc');
    xlabel('Time (samples)');
    legend('Crosscorr', 'Acc data','DetectedSteps')
end
%% insert triggers
%left foot
num=length(AccL.peakSmp);
count=0;
for n = length(EEG.event)+1:length(EEG.event)+num
    count=count+1;
    EEG.event(n).type = 'LeftStep';
    EEG.event(n).latency = AccL.peakSmp(count);
    EEG.event(n).urevent = n;
end

num=length(AccR.peakSmp);
count=0;
for n = length(EEG.event)+1:length(EEG.event)+num
    count=count+1;
    EEG.event(n).type = 'RightStep';
    EEG.event(n).latency = AccR.peakSmp(count);
    EEG.event(n).urevent = n;
end

end



    
            
        

