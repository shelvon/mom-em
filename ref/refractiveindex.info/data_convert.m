delimiter = {'\t',',',' '};

%% for data downloaded from refractiveindex.info
startRow = 1;
formatSpec = '%f%f%[^\n\r]';

%% Open the text file.
% scaleF=1e-9;refName='AlGaAs'
% scaleF=1e-6;refName='AlMcPeak'
scaleF=1e-6;refName='CrJC'


% load n
filename = [refName,'_n.txt'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

wl = dataArray{:, 1}*scaleF;
n = dataArray{:, 2};

% load k
filename = [refName,'_k.txt'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
k = dataArray{:, 2};

% get n complex
nc=n+1i*k;
nr=real(nc);
ni=imag(nc);

%% Write a short table of the exponential function to a text file called AlMcPeak.ref.
format longEng
% dlmwrite([refName,'.ref'],[wl n k], 'delimiter', ' ', 'precision', '%.7e');
dlmwrite([refName,'.ref'],[wl nr ni], 'delimiter', ' ', 'precision', '%.7e');