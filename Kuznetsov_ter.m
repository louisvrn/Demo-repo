clear 
close all
clc

%% Upload data to workspace
% IMPORTANT: must have (x*sensor_sampling_rate) number of samples per
% sensor

fs.NIRS = 50;
fs.EMG = 2048;
load('LT_0.mat');

EMG = load('EMG bip ACC0.mat'); %change numbering for each participant
EMG.raw.data = EMG.EMG_ACC(1,1:535*fs.EMG);
%EMG.raw.data = EMG.EMG_ACC.';

NIRS = load('NIRS0.mat'); %change numbering for each participant
NIRS.raw.data = NIRS.NIRSdata(:,125*fs.NIRS:end);
%NIRS.data.raw = NIRS.data.raw.';

%define time from sampling frequency

NIRS.raw.time = 0:1/fs.NIRS:(length(NIRS.raw.data)-1)/fs.NIRS;
EMG.raw.time = 0:1/fs.EMG:(length(EMG.raw.data)-1)/fs.EMG;

%% Filter EMG signal to get EMG75 

Bpass=designfilt('bandpassiir', 'FilterOrder', 10, 'HalfPowerFrequency1', 100, 'HalfPowerFrequency2', 500, 'SampleRate', fs.EMG);
EMG.filt.data = filter(Bpass, EMG.raw.data);

%% EMG75 'straightening' and averaging over 10s periods

EMG.abs.data = abs(EMG.filt.data);

for i = 1:length(EMG.abs.data)/(fs.EMG*5)
    EMG.avg.data(i) = mean(EMG.abs.data((fs.EMG*(i-1)*5)+1:fs.EMG*i*5));
end

%% Calculate and average HHb

for j = 1:length(NIRS.raw.time)
    i = 0:2;
    HHb.raw.data(j) = sum(NIRS.raw.data(4*i+2,j))/3;
end

for i = 1:length(NIRS.raw.time)
    HHb.adj.data(i) = HHb.raw.data(i)+abs(HHb.raw.data(1));
end

for i = 1:length(NIRS.raw.time)/(fs.NIRS*5)
    HHb.avg.data(i) = mean(HHb.adj.data((fs.NIRS*(i-1)*5)+1:fs.NIRS*(i)*5));
end


%% plotting and comparing EMG and HHb data
% time.avg = 1:length(HHb.avg.data);
% 
% figure()
% subplot(2,1,1)
% plot(time.avg, HHb.avg.data);
% ylabel('HHb')
% subplot(2,1,2)
% plot(time.avg, EMG.avg.data);
% ylabel('EMG')
% xlabel('Time (s)')

%% Calculate ratio HHb/EMG75

ratio.raw.data = HHb.avg.data./EMG.avg.data;
ratio.time = 5:5:5*length(ratio.raw.data);

%% Fit polynomial to approximate data accurately

[ratio.p,~,ratio.mu] = polyfit(ratio.time, ratio.raw.data, 5);
ratio.fit = polyval(ratio.p,ratio.time,[],ratio.mu);

%% Find inflection point 

[data, timestamp] = max(ratio.fit);

LT2 = find(ratio.time>= LT2, 1, 'first');

%% plot final result

figure()
plot(ratio.time, ratio.raw.data);
hold on
plot(ratio.time, ratio.fit, 'r--', ratio.time(timestamp), ratio.fit(timestamp),'g*', ratio.time(LT2), ratio.fit(LT2),'r*', 'MarkerSize', 10)
hold off
xlabel('Time (s)','Interpreter', 'latex')
ylabel('Amplitude (a.u.)','Interpreter', 'latex')
legend({'$\frac{HHb}{EMG_{75}}$', 'Polynomial fit', 'Estimated LT', 'LT'}, 'Interpreter', 'latex')


%% Save results
%change numbering 'Kuznetsov_threshold_1' for each participant in order to
%be able to differentiate participant results
%save('Kuznetsov_timestamp0', 'timestamp')