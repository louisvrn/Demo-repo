clear 
close all
clc

%% upload data to workspace

load('l0.mat'); %change numbering for each participant
lact.time = l0.la_time(2:end); %cut data down to only exercise samples 
lact.data = l0.la(2:end); 

%% find reading errors and approximate blood lactate values at 4 & 6 min

%if lact concentration decreases significantly between two measures,
%identify it as sensor error and re-estimate based on previous and
%following samples
for i = 2:length(lact.data)-1
    if(lact.data(i)<lact.data(i-1)-0.1)
        lact.data(i) = lact.data(i-1)+(lact.data(i+1)-lact.data(i-1))/2;
    end
end

%estimate missing lact concentration values at 4 & 6 min 
lact.four = lact.data(1)+(lact.data(2)-lact.data(1))/2;
lact.six = lact.data(2)+(lact.data(3)-lact.data(2))/2;

%rewrite data and time arrays with all extrapolated values
lact.data = [lact.data(1) lact.four lact.data(2).' lact.six lact.data(3:end).'];
lact.time = [lact.time(1) 4 lact.time(2).' 6 lact.time(3:end).'];

%% fit polynomial to data
%this will allow to have lact estimations per second
p = polyfit(lact.time, lact.data, 5);
fit.time = lact.time(1):1/60:lact.time(end);
fit.data = polyval(p,fit.time);

%% compute thresholds

LT1 = find(fit.data>= 2, 1, 'first')-1;
LT2 = find(fit.data>= 4, 1, 'first')-1;

for i = 1:length(fit.data)-120
    mean1(i) = mean(fit.data(i:i+60));
    mean2(i) = mean(fit.data(i+60:i+120));
    dif(i) = mean1(i)-mean2(i);
end
LTonset = find(dif<-0.5, 1, 'first')-1;

%% Plotting estimation results
figure()
plot(l0.la_time(2:end), l0.la(2:end), 'b.', 'MarkerSize',7, 'linewidth', 4)
hold on
plot(lact.time, lact.data, 'bo');
plot(fit.time, fit.data, 'r--')
plot(fit.time(LTonset), fit.data(LTonset),'mx',fit.time(LT1), fit.data(LT1),'gx', fit.time(LT2), fit.data(LT2),'kx','MarkerSize',10, 'linewidth', 4);
yline(4, '--')
yline(2, '--')
hold off
xlabel('Time (min)','Interpreter', 'latex')
ylabel('Muscle lactate concentration (mM)','Interpreter', 'latex')
legend({'Measured samples','Corrected samples', 'Polynomial fit','Onset in lactate accumulation','LT1', 'LT2'},'Interpreter', 'latex');
yticks(0:2:max(lact.data)+1)

%% Save results
save('LT_0', 'LT1', 'LT2','LTonset')

%% Clear workspace if needed
% clear 
% close all
% clc

%% Load EMG sensor data

datasetsEMG = dir('EMG bip ACC0.mat');

EMG.input.data = load(datasetsEMG.name, '-mat');


%% Constants
% add some constants to the struct const
const.fs = 2048;                          %sample frequency
const.T = 1/const.fs;                     %sample period
raw.data = EMG.input.data.EMG_ACC(1,1:const.fs*540); %raw data
raw.t = (0:length(raw.data)-1)*const.T;   %time of raw signal
load('LT_0.mat')
const.LT2 = LT2/60;                     %LT2 derived from lact data, change for each participant accordingly


%% Filter
% apply a bandpass filter of:
% * order 10
% * low cuttof frequency of 10 Hz
% * high cuttof frequency of 500 Hz
% save filtered data to struct Filtered
Filt.filter = designfilt('bandpassiir', 'FilterOrder', 10, 'HalfPowerFrequency1', 10, 'HalfPowerFrequency2', 500, 'SampleRate', const.fs);
Filt.data = filter(Filt.filter, raw.data);
Filt.time = raw.t;

%% Window split [Razanskas 2015]
% apply the segmentation method used by Razanskas
% save the data to the struct seg

seg.N = 1280;                         %N = estimated window size
i2 = 2;                               %i2 = index for iteration
seg.S = Filt.data;                   %S(t) = filtered data during incremental exercise
seg.t = Filt.time;                   %t = time during incremental exercise
for i = 2:length(seg.S)
    seg.dS(i2) = seg.S(i)-seg.S(i-1); %dS(t), the backward difference of the incremental exercise
    i2 = i2+1;
end

%calculations of V(t), variability function:
i2 = 1;                               %reset index for next iteration
for i = 1:length(seg.dS)
    if i > length(seg.dS)-seg.N
        seg.V(i2) = 0;                                %Prevent Errors at the end of iteration
    else
        seg.V(i2) = sum(abs(seg.dS(i:i+seg.N-1)));    %Equation for V(t)
    end
    i2 = i2+1;
end

%calculations of V_com(t), comparing variability before and after a 
%given sample will allow to determine if it is beginning/end of an a
%activity burst:
i2 = 1;
for i = 1:length(seg.dS)
    if i <= seg.N
        seg.Vcom(i2) = 0;                             %Prevent errors at beginning of iteration
    else
        seg.Vcom(i2) = seg.V(i-seg.N)/seg.V(i);       %Equation for V_com(t)
    end
    i2 = i2+1;
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
end

%% normalized standard deviation per contraction

seg.cycles.stdmax = max(seg.cycles.std);
seg.cycles.nstd = seg.cycles.std/seg.cycles.stdmax;

%% find the first contraction past the LT

seg.lact = find(seg.cycles.start>const.LT2*const.fs*60, 1, 'first');
LT_per_contraction = seg.lact;

%define time axis based on number of contractions
seg.cycles.count = 1:length(seg.cycles.nstd);

%% Extract state specific results (mean & standard deviation of Gaussian curves)

%split and plot observation data into states based on LT
seg.state1.data = seg.cycles.nstd(1:seg.lact);
seg.state2.data = seg.cycles.nstd(seg.lact:end);

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

figure()
seg.state1.prob = histogram(seg.state1.data,'Normalization','probability');
hold on
seg.state2.prob = histogram(seg.state2.data,'Normalization','probability');
hold off
xlabel('STD (per contraction)','Interpreter', 'latex')
ylabel('Probability density','Interpreter', 'latex')
legend({'pre LT', 'post LT'},'Interpreter', 'latex')

figure()
plot(seg.state1.data,y1,'b.', 'MarkerSize', 7);
hold on
plot(seg.state2.data,y2,'r.','MarkerSize', 7);
ylabel('Gaussian fitted distribution','Interpreter', 'latex')
xlabel('Normalized STD (per contraction)','Interpreter', 'latex')
legend({'pre-LT', 'post-LT'},'Interpreter', 'latex')
hold off

nSTD = seg.cycles.nstd;

%% save results
%state-specific mean and sigma values of participant 1 for definition of 
%first iteration of HMM, non state-specific STD contraction values
%for all participants to use as input to train HMM

save('model_parameters_nstd0','nSTD', 'm1', 'm2', 's1', 's2', 'transmat','LT_per_contraction')

%% 
clear
close all
clc

%% Load data to workspace

datasets.parameters = dir('model_parameters_nstd0.mat'); 
for i = 1:length(datasets.parameters)
    data(i) = load(datasets.parameters(i).name, '-mat');
    cycles{i,:} = 1:length(data(i).nSTD);
    nSTD{i}(:,cycles{i,:}) = data(i).nSTD;
end

%% Gaussian fitting
%Now let us introduce the initial Gaussian parameters for Q = 2 hidden states for participant 1 to
%base first iteration of HMM off
%The mean and standard deviation of the Gaussians are used to define
%emission probabilities of the HMM 
Q = 2; 

mu0 = [data(i).m1 data(i).m2];
Sigma0 = [data(1).s1 data(1).s2];

mu0 = reshape(mu0, [1 Q 1]);
Sigma0 = reshape(Sigma0, [1 1 Q 1]);

%Define prior and transition probabilities, the mixmatrix is null
%in this case since we are not using a mixture of Gaussian inputs to 
%represent each state but only one

mixmat0 = [];
prior0 = [1; 0];
transmat0 = data(1).transmat;

%% Iteratively train model
%Finally, let us iteratively improve these parameter estimates using EM
%and the computed STD per contraction

[LL_HMM, prior_HMM, transmat_HMM, mu_HMM, Sigma_HMM, mixmat_HMM] = mhmm_em(nSTD, prior0, transmat0, mu0, (Sigma0).^2, mixmat0, 'max_iter', 20);
 
save('HMM_parameters', 'prior_HMM', 'transmat_HMM', 'mu_HMM', 'Sigma_HMM', 'mixmat_HMM');

%%
% clear 
% close all
% clc
%% Estimating state at every observation for a given data sequence
%load observation data for one participant to test the model
observation = load('model_parameters_nstd0.mat');
data = observation.nSTD;
LT2 = observation.LT_per_contraction; %LT found previously


%%
%load model parameters
model = load('HMM_parameters.mat');
mu = model.mu_HMM;
Sigma = model.Sigma_HMM;
mixmat = model.mixmat_HMM;
transmat = model.transmat_HMM;

%% Use Viterbi algorithm to predict hidden state for every contraction in the input data
B = mixgauss_prob(data, mu, Sigma);
[path] = viterbi_path([1;0], transmat, B);

%% Find contraction at which hidden state transition is predicted to have been made (if at all) 
curr_state = path(end);
transition = find(path>1, 1, 'first');





