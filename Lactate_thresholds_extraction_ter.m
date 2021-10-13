clear 
close all
clc

%% upload data to workspace

load('l0.mat'); %change numbering for each participant
lact.time = l0.la_time(2:end); %cut data down to only efit.timeercise samples 
lact.data = l0.la(2:end);

%% find reading errors and approximate blood lactate values at 4 & 6 min

for i = 2:length(lact.data)-1
    if(lact.data(i)<lact.data(i-1))
        lact.data(i) = lact.data(i-1)+(lact.data(i+1)-lact.data(i-1))/2;
    end
end

lact.four = lact.data(1)+(lact.data(2)-lact.data(1))/2;
lact.fit.time = lact.data(2)+(lact.data(3)-lact.data(2))/2;

lact.data = [lact.data(1) lact.four lact.data(2).' lact.fit.time lact.data(3:end).'];
lact.time = [lact.time(1) 4 lact.time(2).' 6 lact.time(3:end).'];

%% fit polynomial to data
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

%% plot data extrapolation and thresholds
figure()
plot(lact.time, lact.data, 'o');
hold on
plot(fit.time, fit.data, 'r--')
plot(fit.time(LTonset), fit.data(LTonset),'mx',fit.time(LT1), fit.data(LT1),'gx', fit.time(LT2), fit.data(LT2),'kx','MarkerSize',10, 'linewidth', 4);
yline(4, '--')
yline(2, '--')
hold off
xlabel('Time (min)','Interpreter', 'latex')
ylabel('Muscle lactate concentration (mM)','Interpreter', 'latex')
legend({'Samples', 'Polynomial fit','Onset in lactate accumulation', 'LT1', 'LT2'},'Interpreter', 'latex');
yticks(0:2:max(lact.data)+1)

%% Save results
%change numbering 'LT_0' for each participant in order to
%be able to differentiate participant results

save('LT_3', 'LT1', 'LT2','LTonset')






