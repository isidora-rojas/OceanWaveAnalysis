function data = swash_loadTableVData(file,np)
% Reads (layer-averaged) output quantities and store them
% Number of locations is assumed to be known
%
% content of data:
%
% data{1} contains names of variables
% data{2} contains units of variables
% data{3} contains times, if selected
%
% from j = 4 till number of quantities+3
%      data{j,k} contains output variable for each layer k
%
% The format of the TABLE file is fixed (see appendix ??? of SWASH manual)
%
% Author  : Marcel Zijlema
% Date    : September 15, 2010
% Version : 1.0

fid = fopen(file);
if ~(fid>0)
    display([file,'could not be open']);
    data = cell(0,0);
    return
end

for iq=1:5 fgetl(fid); end;                             % skip 5 header lines

[nq] = fscanf(fid,'%i',[1 1]);                          % read number of quantities
fgetl(fid);

for iq=1:nq
   namvar{iq} = fscanf(fid,'%c',[1  8]); fgetl(fid);    % read variable name
   unit{iq}   = fscanf(fid,'%c',[1 16]); fgetl(fid);    % read unit
   kmax(iq)   = fscanf(fid,'%i',[1  1]); fgetl(fid);    % read dimension in vertical
   excv(iq)   = fscanf(fid,'%f',[1  1]); fgetl(fid);    % read execption value
end;

% read rest of table and store data
nt=0;
neof=1;
while neof
   nt=nt+1;
   for ip=1:np
       for iq=1:nq
           vals = fscanf(fid,'%f',[1 kmax(iq)]);
           if isempty(vals), neof=0;, break, end
           vals(vals==excv(iq)) = NaN;
           tabl(iq,ip,nt,1:kmax(iq)) = vals(1:kmax(iq));
       end
       if neof == 0, break, end
   end
end

nt=nt-1;
kmaxm = max(kmax);
data = cell(nq+2,kmaxm);

data{1,1} = namvar;
data{2,1} = unit;

j=0;
for iq=1:nq
    if namvar{iq} == 'Tsec    '
       time = reshape(tabl(iq,1:np,1:nt,1),np,nt);
       data{3,1} = time(1,:)';
    else
       j=j+1;
       for ik=1:kmax(iq)
           dat = reshape(tabl(iq,1:np,1:nt,ik),np,nt);
           data{j+3,ik} = dat';
       end
    end
end

fclose(fid);
