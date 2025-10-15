function [P,F,dof]=crosgk(X,Y,N,M,DT,DW,stats);

% CROSGK   Power cross-spectrum computation, with smoothing in the
%          frequency domain
%
% Usage: [P,F]=CROSGK(X,Y,N,M,DT,DW,stats)
%
% Input:
% X  contains the data of series 1
% Y  contains the data of series 2
% N  is the number of samples per data segment (power of 2)
% M  is the number of frequency bins over which is smoothed (optional),
%    no smoothing for M=1 (default)
% DT is the time step (optional), default DT=1
% DW is the data window type (optional): DW = 1 for Hann window (default)
%                                        DW = 2 for rectangular window
% stats : display resolution, degrees of freedom (optimal, YES=1, NO=0)
%
% Output:
% P contains the (cross-)spectral estimates: column 1 = Pxx, 2 = Pyy, 3 = Pxy
% F contains the frequencies at which P is given
%

%
% Gert Klopman, Delft Hydraulics, 1995
%

if nargin < 4,
  M = 1;
end;

if nargin < 5,
    DT = 1;
end;

if nargin < 6,
  DW = 1;
end;

if nargin < 7,
  stats = 1;
end;

df = 1 / (N * DT) ;

% data window

w = [];
if DW == 1,
  % Hann
  w  = hanning(N);
  dj = N/2;
else,
  % rectangle
  w  = ones(N,1);
  dj = N;
end;
varw = sum (w.^2) / N ;


% summation over segments

nx    = max(size(X));
ny    = max(size(Y));
avgx  = sum(X) / nx;
avgy  = sum(Y) / ny;
px    = zeros(size(w));
py    = zeros(size(w));
Pxx   = zeros(size(w));
Pxy   = zeros(size(w));
Pyy   = zeros(size(w));
ns    = 0;

for j=[1:dj:nx-N+1],

  ns = ns + 1;

  %   compute FFT of signals

  px = X([j:j+N-1]') - avgx;
  px = w .* px ;
  px = fft(px) ;

  py = Y([j:j+N-1]') - avgy;
  py = w .* py ;
  py = fft(py) ;

  % compute periodogram

  Pxx = Pxx + px .* conj(px) ;
  Pyy = Pyy + py .* conj(py) ;
  Pxy = Pxy + py .* conj(px) ;

end;

Pxx = (2 / (ns * (N^2) * varw * df)) * Pxx ;
Pyy = (2 / (ns * (N^2) * varw * df)) * Pyy ;
Pxy = (2 / (ns * (N^2) * varw * df)) * Pxy ;

% smoothing

if M>1,
  w = [];
  w = hamming(M);
  w = w / sum(w);
  w = [w(ceil((M+1)/2):M); zeros(N-M,1); w(1:ceil((M+1)/2)-1)];
  w = fft(w);
  Pxx = fft(Pxx);
  Pyy = fft(Pyy);
  Pxy = fft(Pxy);
  Pxx = ifft(w .* Pxx);
  Pyy = ifft(w .* Pyy);
  Pxy = ifft(w .* Pxy);
end;

Pxx = Pxx(1:N/2);
Pyy = Pyy(1:N/2);
Pxy = Pxy(1:N/2);

% frequencies

F = [];
F = ([1:1:N/2]' - 1) * df;

% signal variance

if DW == 1,
  nn = (ns + 1) * N / 2;
else,
  nn = ns * N;
end;
avgx  = sum (X(1:nn)) / nn;
varx  = sum ((X(1:nn) - avgx).^2) / (nn - 1);
avgy  = sum (Y(1:nn)) / nn;
vary  = sum ((Y(1:nn) - avgy).^2) / (nn - 1);
covxy = sum ((X(1:nn) - avgx) .* (Y(1:nn) - avgy)) / (nn - 1);

m0xx    = (0.5 * Pxx(1) + sum(Pxx(2:N/2-1)) + 0.5 * Pxx(N/2)) * df;
m0yy    = (0.5 * Pyy(1) + sum(Pyy(2:N/2-1)) + 0.5 * Pyy(N/2)) * df;
m0xy    = (0.5 * Pxy(1) + sum(Pxy(2:N/2-1)) + 0.5 * Pxy(N/2)) * df;

%disp(['m0x / varx = ' num2str(m0xx./varx) '  ;  m0y / vary = ' num2str(m0yy./vary) '  ; m0xy / varxy = ' num2str(real(m0xy)./covxy) '  '])


Pxx = Pxx * (varx  / m0xx);
Pyy = Pyy * (vary  / m0yy);
Pxy = Pxy * (covxy / real(m0xy));

P = [Pxx, Pyy, Pxy];

% output spectrum characteristics
dof = floor(2*ns*(M+1)/2/(3-DW));
if stats == 1
fprintf('number of samples used : %8.0f\n', nn);
fprintf('degrees of freedom     : %8.0f\n', floor(2*ns*(M+1)/2/(3-DW)));
fprintf('resolution             : %13.5f\n', (3-DW)*df*(M+1)/2);
end
%
