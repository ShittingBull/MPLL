% ----------------------------------------------------------------------
%   Diese Funktion realisiert ein einfaches AGC-System mit Zeitkonstanten
%   Filtern nach [Zoel95]. 
%
%       agc_test(signal,ta,tr,ampl,thresh,fs)
%
%       signal: Eingangssignal
%       ta: Attackzeit in ms
%       tr: Releasezeit in ms
%       ampl: Zielamplitude (linear)
%       thresh: Schwellwert (dB) 
%       fs: Abtastfrequenz zur Berechnung der Zeitkonstanten 
%        
% H. Riekehof 26.4.2012
% S. Kraft 24.09.2012, kleine Korrekturen
% S. Kraft 18.10.2012, lineare Amplituden in der Schleife
% ----------------------------------------------------------------------

function [y env] = agcfunnew(signal,ta,tr, target_ampl,thresh,fs)

% Zeitkonstantenberechnung nach Zoelzer
Ts = 1/fs;
a = 1 - exp(-2.2*Ts/(ta/1000));
r = 1 - exp(-2.2*Ts/(tr/1000));

thresh = 10^(thresh/20);

% Speicherstellen
peakh = 0;
y = zeros(size(signal));
env = zeros(size(signal));

lastpeak = 0;
env = PDdecoupledsmooth(signal, fs, ta / 1000, tr / 1000);
for n=1:length(signal)
    
    % Spitzenwertgleichrichtung mit Zeitkonstanten
%     s1 = abs(signal(n));
%     if s1 > peakh
%         peakh = (1-a)*peakh + a*s1;
%         lastpeak = n;
%     else        
%         peakh = peakh * ((1-r)^((n-lastpeak)/1000));
%     end
%      
%     % Pegelprozessor
%     env(n) = peakh;
    if env(n) > thresh
        y(n) = (target_ampl / env(n)) * signal(n);    
    else
        y(n) = signal(n);
    end
    
end

end






