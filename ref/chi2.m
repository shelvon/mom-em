clear;close all;
matLabel='GaAs'
% matLabel='AlGaAs'
switch matLabel
	case 'GaAs'
		chi2ZB=[0     0     0   370     0     0
		      	0     0     0     0   370     0
		      	0     0     0     0     0   370]
		
		chi2WZ=[0     0     0     0    42	  0
		      	0     0     0    42     0     0
		       21    21   115     0     0     0]
	case 'AlGaAs'
		chi2ZB=[0     0     0     1     0     0
		      	0     0     0     0     1     0
		      	0     0     0     0     0     1]
end
try
	dlmwrite([matLabel,'ZB','.chi2'],chi2ZB,'delimiter','\t');
    dlmwrite([matLabel,'WZ','.chi2'],chi2WZ,'delimiter','\t');
catch MException
   return;
end 