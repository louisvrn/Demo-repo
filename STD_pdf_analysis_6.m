clear 
close all
clc
%% Load the data
% load the data from the local pc folder
% save the data in the struct raw

datasetsEMG = dir('EMG bip ACC0.mat');

EMG.input.data = load(datasetsEMG.name, '-mat');


%% Constants
% add some constants to the struct const
const.fs=2048;                          %sample frequency
const.T=1/const.fs;                     %sample period
raw.data = EMG.input.data.EMG_ACC(1,1:const.fs*540); %raw data
raw.t=(0:length(raw.data)-1)*const.T;   %time of raw signal
load('LT_0.mat')
const.LT2 = LT2/60;                        %LT2 derived from lact data, change for each participant accordingly


%% Filter
% apply a bandpass filter of:
% * order 10
% * low cuttof frequency of 10 Hz
% * high cuttof frequency of 500 Hz
% save filtered data to struct Filtered
Filtered.filter=designfilt('bandpassiir', 'FilterOrder', 10, 'HalfPowerFrequency1', 10, 'HalfPowerFrequency2', 500, 'SampleRate', const.fs);
Filtered.d=filter(Filtered.filter, raw.data);
Filtered.t=raw.t;

%% Window split [Razanskas 2015]
% apply the segmentation method used by Razanskas
% save the data to the struct seg
seg.N=1280;                         %the estimated window size
i2=2;                               %this is an index for the iteration
seg.S=Filtered.d;%S(t), the filtered data during incremental exercise
seg.t=Filtered.t;              % t, the time during incremental exercise
for i=2:length(seg.S)
    seg.dS(i2)=seg.S(i)-seg.S(i-1); %dS(t), the backward difference of the incremental exercise
    i2=i2+1;
end

%calculations of V(t):
i2=1;                               %reset index for next iteration
for i=1:length(seg.dS)
    if i>length(seg.dS)-seg.N
        seg.V(i2)=0;                                %Prevent Errors at the end of iteration
    else
        seg.V(i2)=sum(abs(seg.dS(i:i+seg.N-1)));    %Equation of V(t)
    end
    i2=i2+1;
end

%calculations of V_com(t):
i2=1;
for i=1:length(seg.dS)
    if i<=seg.N
        seg.Vcom(i2)=0;                             %Prevent errors at beginning of iteration
    else
        seg.Vcom(i2)=seg.V(i-seg.N)/seg.V(i);       %Equation of V_com(t)
    end
    i2=i2+1;
end

%% Split data into contractions using Vcom info

%find all local min and max of Vcom (aka the start and end of every
%contraction)

fig.xlimvar=[300 305];

seg.min = islocalmin(abs(seg.Vcom),'MinProminence',0.75);
seg.max = islocalmax(abs(seg.Vcom),'MinProminence',0.75);

%define and store all starting and ending points of contractions
seg.cycles.start = find(seg.min);
seg.cycles.end = find(seg.max);

%% Plot contraction identification

figure()
plot(seg.t,seg.S,seg.t(seg.min),seg.S(seg.min),'r*',seg.t(seg.max),seg.S(seg.max),'g*')
axis tight
ylabel("Amplitude (a.u.)",'Interpreter', 'latex')
xlabel("Time (s)",'Interpreter', 'latex')
legend({'Filtered EMG', 'Contraction start', 'Contraction end'});
hold on
xlim(fig.xlimvar)
figure()
plot(seg.t,seg.Vcom,seg.t(seg.min),seg.Vcom(seg.min),'r*',seg.t(seg.max),seg.Vcom(seg.max),'g*')
axis tight
ylabel("Vcom",'Interpreter', 'latex')
xlabel("Time (s)",'Interpreter', 'latex')
hold on
xlim(fig.xlimvar)

%% find standard deviation per contraction

%calculate STD filtered EMG amplitude value per identified contraction
for i = 1:length(seg.cycles.start)-1
    seg.cycles.std(i) = std(seg.S(seg.cycles.start(i):seg.cycles.end(i+1)));
    seg.cycles.mean(i) = mean(seg.S(seg.cycles.start(i):seg.cycles.end(i+1)));   
    %seg.Sstd(seg.cycles.start(i):seg.cycles.end(i+1)) = seg.dS(seg.cycles.start(i):seg.cycles.end(i+1));
end

seg.cycles.stdmax = max(seg.cycles.std);
seg.cycles.std = seg.cycles.std/seg.cycles.stdmax;

%find the first contraction past the LT
seg.lact = find(seg.cycles.start>const.LT2*const.fs*60, 1, 'first');
LT_per_contraction = seg.lact;

%define time axis based on number of contractions
seg.cycles.count = 1:length(seg.cycles.std);

%% Exctract state specific results (mean & standard deviation of Gaussian curves)

%split and plot observation data into states based on LT
seg.state1.data = seg.cycles.std(1:seg.lact);
seg.state2.data = seg.cycles.std(seg.lact:end);
% seg.state1.data = seg.Sstd(1:seg.lact*const.fs);
% seg.state2.data = seg.Sstd(seg.lact*const.fs:end);

%Gaussian fitting for both pre and post LT
[m1,s1] = normfit(seg.state1.data);
y1 = normpdf(seg.state1.data,m1,s1);
[m2,s2] = normfit(seg.state2.data);
y2 = normpdf(seg.state2.data,m2,s2);
%% Evaluate transition matrix (only for participant who will be used for first iteration of HMM)

transition11 = 0;
transition22 = 0;
transition12 = 0;
transition21 = 0;

for i = 1:length(seg.cycles.start)-1
    if((seg.cycles.start(i) < seg.lact*const.fs) && (seg.cycles.start(i+1) >= seg.lact*const.fs))
        transition12 = transition12 + 1;
    end
    if((seg.cycles.start(i) >= seg.lact*const.fs) && (seg.cycles.start(i+1) < seg.lact*const.fs))
        transition21 = transition21 + 1;
    end
    if((seg.cycles.start(i) < seg.lact*const.fs) && (seg.cycles.start(i+1) < seg.lact*const.fs))
        transition11 = transition11 + 1;
    end
    if((seg.cycles.start(i) >= seg.lact*const.fs) && (seg.cycles.start(i+1) >= seg.lact*const.fs))
        transition22 = transition22 + 1;
    end
end

transmat = [[transition11 transition12]/(transition11+transition12); [transition21 transition22]/(transition21+transition22)];


%% Plot outputs 

% figure()
% seg.state1.prob = histogram(seg.state1.data,'Normalization','probability');
% hold on
% seg.state2.prob = histogram(seg.state2.data,'Normalization','probability');
% hold off
% xlabel('STD (per contraction)','Interpreter', 'latex')
% ylabel('Probability density','Interpreter', 'latex')
% legend({'pre LT', 'post LT'},'Interpreter', 'latex')
figure()
plot(seg.t,seg.S,seg.t(seg.min),seg.S(seg.min),'r*',seg.t(seg.max),seg.S(seg.max),'g*')
axis tight
ylabel("Amplitude (a.u.)",'Interpreter', 'latex')
xlabel("Time (s)",'Interpreter', 'latex')
legend({'Filtered EMG', 'Contraction start', 'Contraction end'});
hold on
xlim(fig.xlimvar)
hold off
% subplot(2,1,2)
% plot(seg.state1.data,y1,'r.', 'MarkerSize', 7);
% hold on
% plot(seg.state2.data,y2,'b.','MarkerSize', 7);
% ylabel('Gaussian fitted distribution','Interpreter', 'latex')
% xlabel('Normalized STD (per contraction)','Interpreter', 'latex')
% legend({'pre-LT', 'post-LT'},'Interpreter', 'latex')
% hold off

STD = seg.cycles.std;

%% save results
%state-specific mean and sigma values of participant 1 for defintion of 
%first iteration of HMM, non state-specific STD contraction values
%for all participants to use as input to train HMM

save('model_parameters_std0','STD', 'm1', 'm2', 's1', 's2', 'transmat','LT_per_contraction')