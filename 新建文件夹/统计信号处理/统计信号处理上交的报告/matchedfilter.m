function smf = matchedfilter(s,r,fs)

% MATCHEDFILTER calculates the matched filter output
%                 of signal s with replica r.
%
%                 smf = matchedfilter(s,r,fs)
%
%                 input:
%                        s    : the signal
%                        r    : the replica
%                        fs   : the sample frequency
%                 output:
%                        smf  : the matchedfilter output
 
%                 Developed at TNO-FEL
%                 Author: Ren\'e Laterveer,      30/09/1992
%                 some minor changes, RL  7/10/1992
%                 added normalization, RL  26/11/1992


ls = size(s);
if (ls(1) ~= 1),
  s = s';
end
ls = max(ls);

r(r==0) = [];
lr = size(r);
if (lr(1) ~= 1),
  r = r';
end
lr = max(lr);
T = lr/fs;

if (ls < lr),
    error('Replica must have smaller length than signal');
end

%  Normalize replica on energy 1
%r = r/sqrt(r*r');

% Normalize replica on maximum amplitude 1
r = r/abs(max(r));

% determine number of points in fft
p = 2^ceil(log(ls)/log(2));
if p==ls,     %  remove wrap around effects when
  p=2*p;      %  signal is already 2^n in length
end

% correlation by multiplication in frequency domain
y = conj(fft(r,p)).*fft(s,p);

% perform Hilbert transform (stolen from hilbert.m)
% throw away negative frequencies and double positive frequencies
if p ~= 1
%                  freq>0                Nyquest            freq<0
	h = [1; 2*ones(fix((p-1)/2),1); ones(1-rem(p,2),1); zeros(fix((p-1)/2),1)];
	y(:) = y(:).*h;
end

% matched filtered output
smf = abs(ifft(y))';        % envelope
%smf1 = real(ifft(y)');      % real part 
%smf2 = imag(ifft(y)');      % imaginery part 

smf(ls+1:p) = [];           % adjust the output to same length as input signal
smf = smf*2/T/fs;           % normalization according to amplitude signal

%figure;
%plot(smf);

% end
