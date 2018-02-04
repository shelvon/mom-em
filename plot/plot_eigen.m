clear;close all;
%
delimiter = {'\t',' '};
formatSpec = '%f%[^\n\r]';

%
% modelPrefix='../examples/nanodisk/nd'
modelPrefix='../examples/sphere/sphere'
% modelPrefix='../examples/cube/cube'

fileName=[modelPrefix,'-real_eig.txt'];
fileID = fopen(fileName,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'MultipleDelimsAsOne', true,  'ReturnOnError', false);
fclose(fileID);
eig_real = [dataArray{1:end-1}].';

fileName=[modelPrefix,'-imag_eig.txt'];
fileID = fopen(fileName,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'MultipleDelimsAsOne', true,  'ReturnOnError', false);
fclose(fileID);
eig_imag = [dataArray{1:end-1}].';

%%
neig=1:33;
% plot(eig_real([1:neig]),eig_imag([1:neig]),'*');
plot3(eig_real(neig),eig_imag(neig),neig,'-*');grid on;grid minor;
% plot(eig_real,eig_imag,'*');
xlabel('real(\lambda)');ylabel('imag(\lambda)');
eig_imag([1:neig])
eig_real([1:neig])