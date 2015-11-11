function [n_hfo, s_candidates] = find_HFO_peaks(s_data)
%find_HFO_peaks extract high frecuency oscilations (HFO) from any signal.
%Input:
%   -s_data: Any signal in form of vector
%Output:
%   - n_hfo:        It is the number of HFO each s_timeElapsed - Cellarray
%   - s_candidates: Position where starts and ends all HFO detected.
%-------------------------------------------------------------------------
%This function is not validated, please be aware of this.
%-------------------------------------------------------------------------

    global Freq_aux;
    global s_rate;
    global f_signalfraction;
    
    Freq_aux         = [];
    s_rate           = 1024;
    f_signalfraction = 0.85;
    
    s_timeElapsed   = 60;%Time elapsed per sample
    s_time          = 0:1/s_rate:s_timeElapsed;
    s_time(end)     = [];
    s_data          = (s_data -mean(s_data))./norm(s_data);
    s_filter        = designfilt('bandpassiir','FilterOrder',20, ...
                      'PassbandFrequency1',80,'PassbandFrequency2',400, ...
                      'PassbandRipple',0.5, ...
                      'StopbandAttenuation1',60,'StopbandAttenuation2',60, ...
                      'SampleRate', 1024);
                  
    s_candidates    = {};
    n_hfo           = {};
    for i = 1:s_rate*s_timeElapsed:numel(s_data)
        number_hfo_rate = 0;
        fprintf('Time %10d to %10d of %10d... \t', i,...
                    min(i+s_rate*s_timeElapsed-1, numel(s_data)), numel(s_data));
        s_data_temp     = s_data(i:min(i+s_rate*s_timeElapsed-1, numel(s_data)));
        s_data_temp     = s_data_temp./max(s_data_temp);
        s_data_filtered = filter(s_filter, s_data_temp);
        s_data_filtered = s_data_filtered./max(s_data_filtered);
        
        %-----------------Set Standard Desviation
        s_data_std       =  s_data(max(1,(i-s_rate*150)+1):...
                                 min(i+(s_rate*150)-1, numel(s_data)));
        s_time_std       =  0:1/s_rate:numel(s_data_std)/s_rate;
        s_time_std(end)  =  [];
        
        s_std            =  std(abs(hilbert(s_data_filtered)));
                        
        %-----------------Set Threshold
        s_evolvent          = abs(hilbert(s_data_filtered));
        s_threshold         = mean(s_evolvent) + (3*s_std);
        s_data_thresholdRAW = s_evolvent > s_threshold; 
        %7 points are equivalent to 6ms
        
        %-----------------Set points above threshold
        s_uphill_RAW          = find(diff(s_data_thresholdRAW)== 1);
        s_downhill_RAW        = find(diff(s_data_thresholdRAW(...
                                    s_uphill_RAW(1):end))==-1)+...
                                    s_uphill_RAW(1)-1;
        %Just in case starts with downhill
        s_uphill              = s_uphill_RAW;
        s_downhill            = s_downhill_RAW;
        
        %If inter-event-intervale less than 10ms are merge into one event
        for pk = 1:numel(s_downhill_RAW)-1
            if (s_uphill_RAW(pk+1)-s_downhill_RAW(pk))<12   %12 points = 10ms
                s_uphill(s_uphill==s_uphill_RAW(pk+1))      = [];
                s_downhill(s_downhill==s_downhill_RAW(pk))  = [];
            end
        end 
        s_data_thresholdRAW = zeros(size(s_data_thresholdRAW));
        for pk = 1:numel(s_downhill)
            s_data_thresholdRAW(s_uphill(pk):s_downhill(pk)) = 1;
        end
        [s_data_threshold, s_maxLabel] = bwlabel(s_data_thresholdRAW);
        %s_maxLabel shoud be same size as s_downhill and s_uphill
        
        %--------Set possible events greater than 6ms and at least 6 peaks
        %for j = 1:s_maxLabel
        for w =1:numel(s_downhill)    
            %Envelope signal greater than 6ms?
            if sum(s_data_threshold==w) < 7
                s_data_threshold(s_data_threshold==w) = 0;
            else %Yes?
                %Then, inside envelope, there are 6 peaks of filtered signal? 
                s_findpeaks= s_data_filtered(s_uphill(w):s_downhill(w));
                s_sumPeaks = numel(findpeaks(abs(s_findpeaks)));
                v_FT       = [];
                v_diff     = [];
                t_time     = [];
                Freq_aux   = [];
                %t_data    = [];
%                 plot_all(s_time, s_data_temp, s_data_filtered, ...
%                          s_evolvent, s_threshold, s_data_threshold,...
%                          s_data_thresholdRAW, s_findpeaks); 
                if s_sumPeaks < 7
                    s_data_threshold(s_data_threshold==w) = 0;
                else
                    %Then, could be a possible event
                    %to be sure if the event is a HFO, with frecuency-time
                    %representation
                    v_start    = find(diff(s_data_threshold==w)==  1);
                    v_end      = find(diff(s_data_threshold==w)== -1);
                    v_diff     = v_end - v_start +1 ;
                    
                    if (isempty(v_start) && ~isempty(v_end)) %If data start with 1
                        s_data_threshold(s_data_threshold==w) = 0;
                        continue; 
                    end
                    
                    if v_start<800
                        s_data_threshold(s_data_threshold==w) = 0;
                        continue; 
                    end
                    
                    %At this time recording is just starting
                    
                    if ~mod(v_diff,2), v_end=v_end+1; end;
                    %In case number of data is odd
                    if max(v_start-s_rate,1)==1 %v_diff+s_rate*2 must be the entire data
                        t_time  = s_time(1:(s_rate*2)+v_diff);
                        t_data  = s_data_temp(1:(s_rate*2)+v_diff);
                    elseif min(v_end+s_rate, numel(s_time)) == numel(s_time)
                        t_time  = s_time(end-(s_rate*2)+v_diff:end);
                        t_data  = s_data_temp(end-(s_rate*2)+v_diff:end);
                    else
                        t_time  = s_time(max(1,v_start-s_rate):min(v_end+s_rate, end));
                        t_data  = s_data_temp(max(1,v_start-s_rate):min(v_end+s_rate, end));
                    end
                    %Need more points to set frecuency-time rep
                    
                    t_time        = t_time - t_time(floor(numel(t_time)/2)+1);
                    %To set gaussian is needed symetric around zero
                    
                    v_FT          = FrecuencyTime(t_data, t_time);
                    
                    %Set maximal region as a possible event
                    r_thres      = bwlabel(abs(v_FT)>0.6*max(max(abs(v_FT))));
                    r_unique     = unique(r_thres); 
                    r_valueCount = histc(r_thres(:),r_unique); r_valueCount(1)=0;
                    r_regionMax  = find(r_valueCount==max(r_valueCount))-1;
                    r_thres      = (r_thres == r_regionMax(end));%*255; 
                    
                    %Set time fot the image
                    t_time_new = 0:1/s_rate:(numel(t_time)*(1-f_signalfraction)+v_diff)/s_rate;
                    t_time_new = t_time_new + (find(s_data_filtered == ...
                        s_findpeaks(1))-s_rate*(1-f_signalfraction))/s_rate;
                    
                    %Set where starts and ends the possible events on the
                    %filtered signal to see if match with the image
                    p_position      = [find(s_data_filtered == s_findpeaks(1))/1024;...
                                        find(s_data_filtered == s_findpeaks(end))/1024];
                    [~, p_position] = min(abs(repmat(t_time_new,2,1) -...
                                            repmat(p_position,1,numel(t_time_new)))');
                    
                    p_event         = logical(zeros(1, size(r_thres,2)));
                    p_event(p_position(1):min(p_position(2), end)) = 1;
%                     fprintf('size of p_event: %d %d, class: %s\n', size(p_event), class('p_event'));
%                     fprintf('size of max(r_thres): %d %d, class: %s\n', size(max(r_thres)), class('max(r_thres)'));
%                     fprintf('size of p_event.*max(r_thres): %d %d, class: %s\n\n', size(p_event.*max(r_thres)), class('p_event.*max(r_thres)'));

                    p_classification    = max(p_event.*max(r_thres));%Logical
                    %It is one if maximal region match with evolvent on
                    %filtered signal
                    if numel(find(max(r_thres')==1))>70 &&...%all rows greater than
                            numel(find(max(r_thres)==1))<40  %all columns are less than
                        p_classification = 0; %Could be a spike
                    end
                    number_hfo_rate = number_hfo_rate + p_classification;
                    
                    s_data_threshold(s_data_threshold==w) = p_classification;
                    
                    plot_events(s_time, s_data_temp, s_data_filtered, ...
                         s_evolvent, s_threshold, s_data_threshold,...
                         s_findpeaks, v_FT, t_time_new, r_thres, p_classification);            
                     
                end
%                 plot_all(s_time, s_data_temp, s_data_filtered, ...
%                          s_evolvent, s_threshold, s_data_threshold,...
%                          s_data_thresholdRAW, s_findpeaks, v_FT, t_time, v_diff);
            end
        end
        fprintf('Possibles HFO found=%d\n', number_hfo_rate);
        n_hfo{end+1}           = number_hfo_rate;
        
        if isempty(find(diff(s_data_threshold) ==  1)), continue; end;
        
        s_candidates{end+1, 1} = find(diff(s_data_threshold) ==  1);
        s_candidates{end, 2}   = find(diff(s_data_threshold) == -1);
    end
    s_candidates    = sort(cell2mat(s_candidates')');  
    
end



function plot_all(s_time, s_data_temp, s_data_filtered, s_evolvent,...
                  s_threshold, s_data_threshold, s_data_thresholdRAW,...
                  s_findpeaks, v_FT, t_time, v_diff)
    global Freq_aux;
    global s_rate;
    global f_signalfraction;
    
    ax1 = subplot(611);  
        plot(s_time, s_data_temp, 'k');ylim([-1 1]);hold on;
        line(repmat(find(s_data_filtered == s_findpeaks(1))/1024,1,2),... 
            [-1 1], 'LineStyle',':', 'Color',[1 0 1], 'LineWidth', 2);
        line(repmat(find(s_data_filtered == s_findpeaks(end))/1024,1,2),... 
            [-1 1], 'LineStyle',':', 'Color',[1 0 1], 'LineWidth', 2);
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_temp)*0.9,'A',...
            'Color', 'red', 'FontSize',15); 
        hold off;
                         
    ax2 = subplot(612);  
        plot(s_time, s_data_filtered, 'k');hold on;
        plot(s_time, s_evolvent, 'r');ylim([-1 1]);
        plot(s_time, ones(1,numel(s_time))*s_threshold, '--b');
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_filtered)*0.9,'B',...
            'Color', 'red', 'FontSize',15); 
        hold off;
                         
    ax3 = subplot(613);  
        plot(s_time, s_data_temp, 'k');hold on;ylim([-1 1]);
        plot(s_time, s_data_thresholdRAW, '-r'); hold off;
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_temp)*0.9,'C',...
            'Color', 'red', 'FontSize',15); 
                         
    ax4 = subplot(614);  
        plot(s_time, s_data_temp, 'k');hold on;ylim([-1 1]);
        plot(s_time, s_data_threshold, '-r'); hold off;
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_temp)*0.9,'D',...
            'Color', 'red', 'FontSize',15);
        if ~isempty(Freq_aux)
            line(repmat((find(s_data_filtered == ...
                s_findpeaks(1))-s_rate*(1-f_signalfraction))/s_rate,1,2),... 
                [-1 1], 'LineStyle',':', 'Color',[1 1 0], 'LineWidth', 2);
            line(repmat((find(s_data_filtered == ...
                s_findpeaks(end))+s_rate*(1-f_signalfraction))/s_rate,1,2),... 
                [-1 1], 'LineStyle',':', 'Color',[1 1 0], 'LineWidth', 2);                                                          
        end
                         
    linkaxes([ax1,ax2,ax3,ax4], 'x');
    ax1.XLim = [s_time(find(s_data_filtered == s_findpeaks(1))) - ...
                (s_rate*(1-f_signalfraction+0.05))/1024,...
                 s_time(find(s_data_filtered == s_findpeaks(end)))...
                 +(s_rate*(1-f_signalfraction+0.05))/1024];
          
    subplot(615);  
        findpeaks(abs(s_findpeaks),s_rate);
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_temp)*0.9,'E',...
            'Color', 'red', 'FontSize',15);
    
    subplot(616);  
        freezeColors;
        if isempty(Freq_aux)
            image(zeros(100,100));
            text(47,50,'X', 'Color', 'red', 'FontSize',40);
            set(gca,'xtick',[],'ytick',[],'layer','bottom','box','on')
            colormap(jet);xlabel('Time [s]');ylabel('Freq [Hz]');
        else
            t_time_new = 0:1/s_rate:(numel(t_time)*(1-f_signalfraction)+v_diff)/s_rate;
            t_time_new = t_time_new + (find(s_data_filtered == ...
                s_findpeaks(1))-s_rate*(1-f_signalfraction))/s_rate;
            imagesc(t_time_new, Freq_aux, abs(v_FT),...
                [min(min(abs(v_FT))) 0.6*max(max(abs(v_FT)))]);
            text(max(s_time)*0.95, max(Freq_aux)*0.92,...
                'F', 'Color', 'red', 'FontSize',15);
            colormap(jet);xlabel('Time [s]');ylabel('Freq [Hz]');
        end
    drawnow;
end

function plot_events(s_time, s_data_temp, s_data_filtered, s_evolvent,...
                  s_threshold, s_data_threshold, s_findpeaks, v_FT,...
                  t_time_new, r_thres, p_classification)
    global Freq_aux;
    global s_rate;
    global f_signalfraction;
    
    ax1 = subplot(511);  
        plot(s_time, s_data_temp, 'k');ylim([-1 1]);hold on;
        line(repmat(find(s_data_filtered == s_findpeaks(1))/1024,1,2),... 
                [-1 1], 'LineStyle',':', 'Color',[1 0 1], 'LineWidth', 2);
        line(repmat(find(s_data_filtered == s_findpeaks(end))/1024,1,2),... 
                [-1 1], 'LineStyle',':', 'Color',[1 0 1], 'LineWidth', 2);
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_temp)*0.9,'A',...
               'Color', 'red', 'FontSize',15); 
        hold off;
                         
    ax2 = subplot(512);  
        plot(s_time, s_data_filtered, 'k');hold on;
        plot(s_time, s_evolvent, 'r');ylim([-1 1]);
        plot(s_time, ones(1,numel(s_time))*s_threshold, '--b');
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_filtered)*0.9,'B',...
            'Color', 'red', 'FontSize',15); 
        hold off;
                                            
    ax3 = subplot(513);  
        plot(s_time, s_data_temp, 'k');hold on;ylim([-1 1]);
        plot(s_time, s_data_threshold, '-r'); hold off;
        ylabel('EEG [uV]');
        text(max(s_time)*0.95,max(s_data_temp)*0.9,'C',...
            'Color', 'red', 'FontSize',15);
        if ~isempty(Freq_aux)
            line(repmat((find(s_data_filtered == ...
                s_findpeaks(1))-s_rate*(1-f_signalfraction))/s_rate,1,2),... 
                [-1 1], 'LineStyle',':', 'Color',[1 1 0], 'LineWidth', 2);
            line(repmat((find(s_data_filtered == ...
                s_findpeaks(end))+s_rate*(1-f_signalfraction))/s_rate,1,2),... 
                [-1 1], 'LineStyle',':', 'Color',[1 1 0], 'LineWidth', 2);                                                          
        end
                         
    linkaxes([ax1,ax2,ax3], 'x');
    ax1.XLim = [s_time(find(s_data_filtered == s_findpeaks(1)))- ...
            (s_rate*(1-f_signalfraction+0.05))/1024,...
            s_time(find(s_data_filtered == s_findpeaks(end)))+...
            (s_rate*(1-f_signalfraction+0.05))/1024];
                  
    subplot(514);  
        freezeColors;
        if isempty(Freq_aux)
            image(zeros(100,100));
            text(47,50,'X', 'Color', 'red', 'FontSize',40);
            set(gca,'xtick',[],'ytick',[],'layer','bottom','box','on')
            colormap(jet);ylabel('Freq [Hz]');
        else
            imagesc(t_time_new, Freq_aux, abs(v_FT),...
                [min(min(abs(v_FT))) 0.8*max(max(abs(v_FT)))]);
            
            text(max(s_time)*0.95, max(Freq_aux)*0.92,...
                'D', 'Color', 'red', 'FontSize',15);
            colormap(jet);ylabel('Freq [Hz]');
        end
    
    subplot(5,1,5); 
        if isempty(Freq_aux)
            image(zeros(100,100));
            text(47,50,'X', 'Color', 'red', 'FontSize',40);
            set(gca,'xtick',[],'ytick',[],'layer','bottom','box','on')
            colormap(jet);xlabel('Time [s]');ylabel('Freq [Hz]');
        else
            imagesc(t_time_new, Freq_aux,cat(3,r_thres, r_thres, r_thres*0));
            hold on;
            text(max(s_time)*0.95, max(Freq_aux)*0.92,...
                'E', 'Color', 'red', 'FontSize',15);
            xlabel('Time [s]');ylabel('Freq [Hz]');
            if p_classification
                text(min(t_time_new)+(max(t_time_new)-min(t_time_new))*0.95,...
                    max(Freq_aux)*0.8, '✓', 'Color', 'green', 'FontSize',50);
            else
                text(min(t_time_new)+(max(t_time_new)-min(t_time_new))*0.95,...
                    max(Freq_aux)*0.8, 'X', 'Color', 'red', 'FontSize',40);
            end
            hold off;
        end
end

function channel = FrecuencyTime(Signal, Time)
%Vector2FrecuencyTime returns a signal with frecuency time domain (array).
%Requires a vector of signal and a vector of time.

    global Freq_aux
    global s_rate
    global f_signalfraction;
    
    minF        = 80;
    maxF        = 400;
    Num_cicles  = 6;
    resF        = 1;

    Freq_aux = minF:resF:maxF;

    Signal      = Signal - mean(Signal);
    Time_aux    = linspace(min(Time),max(Time),numel(Signal));
    N           = numel(Time_aux);
    channel     = zeros(numel(Freq_aux), N);
    Signal_fft  = fft(Signal);
    
    %Time_aux    = Time_aux + max(Time_aux);
    %Uncomment if validation is enabled
    
    for iter=1:numel(Freq_aux)
        cycle_s = 1/Freq_aux(iter);
        std_s   = cycle_s * Num_cicles/2;
        win1    = exp(-0.5*Time_aux.*Time_aux/(std_s.*std_s));
        win     = win1.*exp(1i*2*pi*Freq_aux(iter)*Time_aux);
        win     = win./norm(win);
        winFFT  = fft(win);
%        fprintf('Applying Frecuency-Time domain for F=%d, extracting evolve signal\n', Freq_aux(iter));
        channel(iter,:) = abs(ifft(winFFT.*Signal_fft));
        warning('off', 'all');
%        channel(iter,:) = env_secant(Time_aux, channel(iter,:), 25,'top');
%         ax1=subplot(7,1,1);plot(Time_aux,win1);title('Gaussian');
%         ax2=subplot(7,1,2);plot(Time_aux,abs(win));title('Gaussian by Cos');
%         ax7=subplot(7,1,3);plot(Time_aux,Signal);title('Original');
%         ax3=subplot(7,1,4);plot(Time_aux,abs(Signal_fft)); title('FFT original');
%         ax4=subplot(7,1,5);plot(Time_aux,abs(winFFT)); title('FFT gauss by cos, movil');
%         ax5=subplot(7,1,6);plot(Time_aux,(abs(winFFT).*abs(Signal_fft))); title('Dot product last two');
%         ax6=subplot(7,1,7);plot(Time_aux,10*abs(ifft(winFFT.*Signal_fft))); title(sprintf('Frecuency = %d',Freq_aux(iter))); 
%         linkaxes([ax1, ax2, ax3, ax4, ax5, ax6, ax7], 'x');
%         pause;
%         drawnow;
%         break;
    end
    channel = channel(:,(s_rate*f_signalfraction):end-(s_rate*f_signalfraction));%Take off data due to edge trouble
    channel = channel./max(channel(:));
    %figure;subplot(2,1,1);plot(Time_aux,Signal);subplot(2,1,2);imagesc(abs(channel));colormap(jet);colorbar;
end