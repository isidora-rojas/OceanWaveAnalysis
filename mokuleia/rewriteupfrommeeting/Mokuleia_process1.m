addpath 'codes'
addpath '/Volumes/LaCie/MARKMOKU91511/MOKUPROCESS'

%A1  Seabird
load MoA18411

%convert pressure from psi to meters of water
%Assume constant water temp of 25°C for now for density calculation
%Redo later with measured temp if available 
%Assume constant atmospheric pressure of 14.7 psi for now
%Redo later with observed atmospheric pressure 

pm = psi2m(pclip,25,14.7);
clear pclip
%43180 seconds of data are measured every 12 hours
%20 second gap to download data to memory
%insert NaN's here for the gap and turn into one
%continuous time series
pp = nan(43200,184);
pp(1:43180,:) = pm;
clear pm
pp = pp(:);
tt = (0:length(pp)-1)'/60/60/24+tclip(1);
clear tclip
pp = pp(1:end-20);
tt = tt(1:end-20);



