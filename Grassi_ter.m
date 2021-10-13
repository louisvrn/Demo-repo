clear  
close all
clc

%% Upload data to workspace
%use only exercise data

fs.NIRS = 50;
struct = load('NIRS0.mat'); %change numbering for each participant
data.raw = struct.NIRSdata(:,fs.NIRS*60:end);
%data.raw = data.raw.';
time.raw = 0:1/fs.NIRS:(length(data.raw)-1)/fs.NIRS;


%% Computing features

% Set all values as offset from first data point for each set of
% measurements

for i = 1:12
    for j = 1:length(data.raw)
        data.adj(i,j) = data.raw(i,j) + abs(data.raw(i,1));
    end
end

for j = 1:length(data.adj)
    i = 0:2;
    O2Hb.raw.data(j) = sum(data.adj(4*i+1,j))/3;
    HHb.raw.data(j) = sum(data.adj(4*i+2,j))/3;
end

%compute and filter HbDiff data
HbDiff.raw.data = O2Hb.raw.data - HHb.raw.data;
HbDiff.filt.data = movmean(HbDiff.raw.data, 5*fs.NIRS);


%% Plotting

% figure()
% plot(time.raw, HbDiff.filt.data)
% xlabel('Time (s)')
% ylabel('HbDiff (a.u)')

%% Take avg data points

%taking average HbDiff over 1 min during exercise

for i = 1:length(time.raw)/(fs.NIRS*60)
    HbDiff.min.data(i) = mean(HbDiff.filt.data((fs.NIRS*(i-1)*60)+1:(fs.NIRS*60*i)));
end

time.min = 1:length(time.raw)/(fs.NIRS*60);

%sample HbDiff filtered every second

% for i = 0:(length(time.raw)/(fs.NIRS))-1
%     HbDiff.samp.data(i+1) = HbDiff.filt.data(fs.NIRS*i+1);
% end
% 
% time.sec = 1:length(time.raw)/(fs.NIRS);

%% polynomial fitting without averaging

[p.raw,~,mu.raw] = polyfit(time.raw, HbDiff.filt.data, 5);
fit.raw = polyval(p.raw, time.raw, [], mu.raw);

p.deriv = polyder(p.raw);
p.deriv2 = polyder(p.deriv);
fit.deriv2 = polyval(p.deriv2, time.raw, [], mu.raw);

for i = 1:length(fit.deriv2)-1
    if((fit.deriv2(i) > 0) && (fit.deriv2(i+1) < 0))
        threshold.raw(i) = 1;
    end
    if((fit.deriv2(i) < 0) && (fit.deriv2(i+1) > 0))
        threshold.raw(i) = 1;
    end
end
  
threshold.raw = (find(threshold.raw, 1, 'first'));

figure()
subplot(2,1,1)
plot(time.raw, HbDiff.filt.data);
hold on
plot(time.raw, fit.raw, 'r--')
plot(time.raw(threshold.raw), HbDiff.filt.data(threshold.raw), 'kx','MarkerSize',10, 'linewidth', 4)
hold off
xlabel('Time (s)','Interpreter', 'latex')
ylabel('Amplitude (a.u.)','Interpreter', 'latex')
legend({'Sampled HbDiff', 'Polynomial fit','Inflection point'},'Interpreter', 'latex');

subplot(2,1,2)
plot(time.raw, fit.deriv2, 'm--');

threshold.raw = (threshold.raw/fs.NIRS)+60;

%% Polynomial fitting with 1 min averaging

[p.min,~,mu.min] = polyfit(time.min, HbDiff.min.data, 5);
fit.min = polyval(p.min, time.min, [], mu.min);

p.deriv = polyder(p.min);
p.deriv2 = polyder(p.deriv);
fit.deriv2 = polyval(p.deriv2, time.min, [], mu.min);

for i = 1:length(fit.deriv2)-1
    if((fit.deriv2(i) > 0) && (fit.deriv2(i+1) < 0))
        threshold.min(i) = 1;
    end
    if((fit.deriv2(i) < 0) && (fit.deriv2(i+1) > 0))
        threshold.min(i) = 1;
    end
end
  
threshold.min = (find(threshold.min, 1, 'first'));

figure()
subplot(2,1,1)
plot(time.min, HbDiff.min.data, 'o');
hold on
plot(time.min, fit.min, 'r--')
plot(time.min(threshold.min), HbDiff.min.data(threshold.min), 'kx','MarkerSize',10, 'linewidth', 4)
hold off
xlabel('Time (s)','Interpreter', 'latex')
ylabel('Amplitude (a.u.)','Interpreter', 'latex')
legend({'Averaged HbDiff', 'Polynomial fit','Inflection point'},'Interpreter', 'latex')

subplot(2,1,2)
plot(time.min, fit.deriv2, 'm--');

threshold.min = (threshold.min*60)+60;


%% Save results
%change numbering 'Grassi_timestamps_0' for each participant in order to
%be able to differentiate participant results
save('Grassi_timestamps0','threshold');




