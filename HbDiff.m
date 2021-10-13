function [rawHbDiff, filtHbDiff] = HbDiff(data, freq, avg)
validateattributes(data,{'double'}, {'size', [14, 30001]})
%Computing HbDiff from given data
%   First the oxy and deoxygenated concentrations of hemoglobine are
%   calculated based on the data received by the three receivers (stored in
%   'data'). The sampling freq 'freq' and averaging time period 'avg' are
%   also required inputs.
%   The Diff between these two values is then calculated to find HbDiff,
%   which is then filtered to eliminate too strong fluctuations. 

for i = 1:12
    for j = 1:length(data)
        data(i,j) = data(i,j) + abs(data(i,1));
    end
end

for j = 1:length(data)
    i = 0:2;
    rawO2Hb(j) = sum(data(4*i+1,j))/3;
    rawHHb(j) = sum(data(4*i+2,j))/3;
end

%compute and filter HbDiff data
rawHbDiff = rawO2Hb - rawHHb;
filtHbDiff = movmean(rawHbDiff, 5*2048);

end

