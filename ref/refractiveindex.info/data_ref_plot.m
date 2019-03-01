% filename = 'AlMcPeak.ref';
% filename = 'AlGaAs.ref';
filename = 'agjc.ref';
% filename = 'CrJC.ref';
delimiter = {'\t',',',' '};
startRow = 1;

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
ref = [dataArray{1:end-1}];

%% plot the data
figure;hold on;
plot(ref(:,1),ref(:,2),'.')
plot(ref(:,1),ref(:,3),'.')