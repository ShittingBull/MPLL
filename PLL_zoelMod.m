
% Standard  ZPLL core as described in 
%
% PLL-BASED PITCH DETECTION AND TRACKING FOR AUDIO SIGNALS, 
% Zölzer, Sankarababu, Möller (2012). 

function [F0, f_osc, yc, ys, xd, xd_lp] = PLL_zoelStd(samples,fs,fc_Pll, Kd, fCenter)

    % init variables with zeros
        nsamples = length(samples);

        F0(nsamples)    = 0;
        yc(nsamples)    = 0;
        ys(nsamples)    = 0;
        f_osc(nsamples) = 0;
        xd_lp(nsamples) = 0;
        xd(nsamples)    = 0;
        phasecounter = 0;
        filtPll(5)     = 0;
        fRangeOsc      = fCenter;
        KO = 2;
        Kd
        %Kd = 2000;
    % init filter coefficients
    %loop filter 
        [B_Pll,A_Pll] = getFilterCoeff(fc_Pll,fs,1/3,'lowpass');

    % init oscillators
        osc_state1(1) = 1;
        osc_state1(2) = 0; 
        fCenter
        %[b,a] = ellip(1,3,30,[(fCenter - pllIntervall) / (fs/2) ,  (fCenter + pllIntervall) / (fs/2)]);
        %samples = filter (b,a,samples);                                 
                                          
    cutoffs(1) = 0;
    for i = 2:nsamples
        x_in = samples(i);
       
        %xd(i) = x_in * yc(i-1) * Kd;
        xd(i) = x_in * yc(i-1) * Kd ;

        % LP filter in Pll
        [ xd_lp(i), filtPll ] = applyFilter( xd(i), B_Pll, A_Pll, filtPll );
        
        
        % get F0
        %F0(i) = abs(2*xd_lp(i) + 1);
        F0(i) = fRangeOsc + KO * xd_lp(i);
%         if i>1
%             F0(i) =  F0_(i) * (1 - 0.98) + 0.98 * F0(i-1);
%         else
%             F0(i) = F0_(i) * (1 - 0.98);
%         end
        % generate oscillator control signal
        % alpha weighs the direct path signal in relation to the LP
        % signal
        alpha = 0;
        %alpha = 0.35;
        f_osc(i) = KO * (xd_lp(i)*(1-alpha) + xd(i)*alpha) + fRangeOsc;       
        %f_osc(i) = (xd_lp(i)*(1-alpha) + xd(i)*alpha)*2 + 1;

        % PLL oscillator
        %oscillator states holds former outputsamples state(1) real
        %state(2) complex
        phasecounter = mod((phasecounter + f_osc(i)/fs * 2 * pi), 2*pi);
        yc(i) = cos(phasecounter + pi / 2);
        %modified for lock indicator --> gardner 
        %[yc(i), ys(i), osc_state1] = NCO_complex(f_osc(i) / fs, osc_state1);        
        %yc(i) = yc(i) *  x_in;
    end
% 
%     figure;
%     subplot 211
%     BS = 8192;
%     plotSpectrogram(samples, BS, 512, fs, 'ylim', [1 600]);
%     hold on;
%     plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end), 'k', 'LineWidth',2);
%     title('Standard ZPLL');   
    %figure(5);
    %plotSpectrogram(out, 4096, 512, fs, 'ylim', [1 1000]);
end