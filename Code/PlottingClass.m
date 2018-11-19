classdef PlottingClass
    % a class def file which has to be the same name as the file
   properties
      % the path to the directory
      path
      % the filename of the .mat file
      filename
      % the entire filepath
      filepath
      % the data points of the .mat
      rawdata
      % name for saving plots
      newname
      % location for the images
      imagedir
      % the range of the data to view
      Range
      % the locations of the red vertical lines
      left
      right
      % for whether or not to open figures
      figs
      % the data of the microphones
      R1
      R2
      R3
   end
   
   methods
       %% Range Cook and Weber plots
       function setzone = setting(obj, data)
           % create a time array
           time = linspace(0, 20, 1200000);
           % check if the starting range specified
           if obj.Range(1) ~= 0
               start = obj.Range(1) * length(time);
           else 
               start = 1;
           end
           % check if ending range specified
           if obj.Range(2) ~= 1
               stop = obj.Range(2) * length(time);
           end
           % set the range of the input data point, eg. R1, R2 or R3
           if exist('stop') == 1
               zz = data(start:stop);
           else
               zz = data(start:end);
           end
           % return variable with range to plot
           setzone = zz;
       end
       %% Plot stemplots
       function [] = stemplots(~, freq, ampspec, color, xtop, ytop, thistitle, thisxlabel, thisylabel)
           stem(freq,ampspec, color);
           xlim([0 xtop]);
           ylim([0 ytop]);
           title(thistitle,'fontsize',15);
           xlabel(thisxlabel,'fontsize',15);
           ylabel(thisylabel,'fontsize',15);
           %set(gca,'Ydir','reverse')
           set(gca,'FontSize',12);
       end
       
       %% Plot db plots
       function [] = dbplots(~, freq, db_spec, color, ybottom, xtop, ...
               thistitle, thisxlabel, thisylabel)
           plot(freq, db_spec, color);
           ylim([ybottom 0]);
           xlim([0 xtop]);
           title(thistitle,'fontsize',15);
           xlabel(thisxlabel,'fontsize',15);
           ylabel(thisylabel,'fontsize',15);
           set(gca,'FontSize',12);
       end
 
       function [freq, amp_spec, db_spec, n2] = fourier(~, variable, FoldingFreq)
           Yy = fft(variable);
           n2 = length(variable)/2;
           n2 = ceil(n2);
           freq = linspace(0, FoldingFreq, n2);
           amp_spec = abs(Yy)/n2;
           db_spec = mag2db(amp_spec);
       end
       %% Fourier Domain - Cook and Weber 
       function [] = Perform(obj)     
           L = length(obj.R1);
           % the length will always have a certain delta t associated
           timestep = 1200 / L;
           
           % sample rate 
           fs = 1 / timestep;        % Hz !
           FoldingFreq = fs / 2;
           
           thisone = figure;
           set(thisone, 'Visible', obj.figs);
           
          
           
           %% Linear Magnitude (C&W)
           
           subplot(2,2,1)
           % pass the object into setting method
           variable = obj.setting(obj.R1);
           
           % pass into fourier calculator method
           [freq1, amp_specx, db_spec1, n1] = ...
               obj.fourier(variable, FoldingFreq);
           
           % plot the stemplots
           thistitle = 'Microphone 1: Roof';
           obj.stemplots(freq1, amp_specx(1:n1), '-b', 100, 0.05, ...
               thistitle , 'Frequency (Hz)', 'Linear Magnitude')
      
           subplot(2,2,2)
           % pass the object into setting method
           variable2 = obj.setting(obj.R2);
           
           % pass into fourier calculator method
           [freq2, amp_specy, db_spec2, n2] = ...
               obj.fourier(variable2, FoldingFreq);
           
           % plot the stemplot
           thistitle = 'Microphone 2: South';
           obj.stemplots(freq2, amp_specy(1:n2), '-g', 100, 0.05, ...
               thistitle , 'Frequency (Hz)', 'Linear Magnitude')

           subplot(2,2,[3, 4])
           % pass the object into setting method
           variable3 = obj.setting(obj.R3);
           
            % pass into fourier calculator method
           [freq3, amp_specz, db_spec3, n3] = ...
               obj.fourier(variable3, FoldingFreq);
           
          % plot the stemplot
           thistitle = 'Microphone 3: North';
           obj.stemplots(freq3, amp_specz(1:n3), '-r', 100, 0.05, ...
               thistitle, 'Frequency (Hz)', 'Linear Magnitude')
           
           % save the figure
           thisplot = ' - freqdomain';
           obj.Putintofiles(thisone, thisplot);          
           %% Decibal Magnitude (C&W)
           thisone2 = figure;
           set(thisone2, 'Visible', obj.figs);
           
           
           subplot(2,2,1)
           % pass data into plotdb
           thistitle = 'Microphone 1: Roof';
           obj.dbplots(freq1, db_spec1(1:n1), '-b', -200, 100, ...
               thistitle, 'Frequency (Hz)', 'Decibal Magnitude (dB)')
      
           subplot(2,2,2)
           % pass data into plotdb
           thistitle = 'Microphone 2: South';
           obj.dbplots(freq2, db_spec2(1:n2), '-g', -200, 100, ...
               thistitle, 'Frequency (Hz)', 'Decibal Magnitude (dB)')

           subplot(2,2,[3, 4])
           % pass data into plotdb
           thistitle = 'Microphone 3: North';
           obj.dbplots(freq3, db_spec3(1:n3), '-r', -200, 100, ...
               thistitle, 'Frequency (Hz)', 'Decibal Magnitude (dB)')
           
           thisplot = ' - db';
           obj.Putintofiles(thisone2, thisplot);         
       end      
      %% Lowpass filter (Butterworth) - Elbing 
      function [] = filtering(obj, s)
          
          fs = 1000;            % sample frequency
          
          dT = 30;                %segment period in seconds
          Ov = 0.10;              %overlap in percentage of segment period
          ONseg = 0;              %1=plot every segment; 0=only plot average for each mic

          Sfilter = 1;

         if Sfilter == 1
             fcut = 200;
             order = 30;
         end
          
         h = figure;
         set(h, 'Visible', obj.figs);

         a=floor(obj.Range(1)*(size(s,1)-1)/(1-0)+1);
         b=ceil(obj.Range(2)*(size(s,1)-1)/(1-0)+1);
         s = s(a:b,:);
         t=(0:1/fs:size(s,1)/fs-1/fs)'; 
         
                 % Lowpass filter (Butterworth)
        if Sfilter == 0
            sf = s;
        end
        if Sfilter == 1
            [b,a] = butter(order,fcut(1)/(fs/2),'low');
            sf = filter(b,a,s);
        end

        for mics = 1:size(s,2)
            for seg = 1:floor(size(s,1)/(dT*fs-Ov*dT*fs))
                 aa = (seg-1)*dT*fs*(1-Ov)+1; bb = aa+dT*fs-1;
                Sf =sf(aa:bb,mics);

                % power spectrum calculation
                N = length(Sf);
                xdft = fft(Sf);
                xdft = xdft(1:N/2+1);
                psdx = (1/(fs*N)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                freq = 0:fs/length(Sf):fs/2;

                % save the power spectrum for each 
                if seg == 1
                    RR(:,1) = freq;
                end
                RR(:,seg+1) = psdx; 
       
            % plot individual spectra
            if mics == 1
                mkr = '--m';
            end
            if mics == 2
                mkr = '--y';
            end
            if mics == 3
                mkr = '--c';
            end

            if ONseg == 1
                semilogx(freq,10*log10(psdx),mkr)
                %         semilogx(freq,20*log10(psdx/10^-6),mkr)
                grid on
                hold on
            end
        end
        if mics == 1
            mkr = '-b';
        end
        if mics == 2
            mkr = '-g';
        end
        if mics == 3
            mkr = '-r';
        end

        semilogx(RR(:,1),10*log10(mean(RR(:,2:size(RR,2)),2)),mkr,'LineWidth',2)
        %     semilogx(RR(:,1),20*log10(mean(RR(:,2:size(RR,2)),2)/10^-6),'-r','LineWidth',2)
        hold on
        clearvars RR
    end

    xlabel('Frequency (Hz)','fontsize',15)
    ylabel('Power/Frequency (dB/Hz)','fontsize',15)
    title('Spectra', 'fontsize',15)

    plot(obj.right, [-120,0],'--r')
    plot(obj.left, [-120,0],'--r', 'HandleVisibility','off')
    legend({'Microphone 1: Roof', 'Microphone 2: South', 'Microphone 3: North', 'Range'}, 'location', 'northwest')
    set(gca,'FontSize',12);

    axis([.5 100 -120 0])

    
    thisplot = ' - spectra';
    obj.Putintofiles(h, thisplot);
      end
      %% Raw Data - Cook and Hartzler
      function [] = Data(obj)
        time = .001:.001:1200;

        l = figure;
        set(l, 'Visible', obj.figs);

        subplot(2,2,1);
        plot(time, obj.R1, '-b')
        title('Microphone 1: Roof', 'fontsize',15)
        ylim([-1 1])
        xlabel('Time (seconds)','fontsize',15)
        ylabel('Pressure [Pa]','fontsize',15)
        set(gca,'FontSize',12);

        subplot(2,2,2);
        plot(time, obj.R2, '-g')
        title('Microphone 2: South', 'fontsize',15)
        ylim([-1 1])
        xlabel('Time (seconds)','fontsize',15)
        ylabel('Pressure [Pa]','fontsize',15)
        set(gca,'FontSize',12);

        subplot(2,2,[3,4]);
        plot(time, obj.R3, '-r')
        title('Microphone 3: North', 'fontsize',15)
        ylim([-1 1])
        xlabel('Time (seconds)','fontsize',15)
        ylabel('Pressure [Pa]','fontsize',15)
        set(gca,'FontSize',12);

        thisplot = 'rawdata - ';
        obj.Putintofiles(l, thisplot);
        % end the plotting function
      end
      
      %% Save the figures
      function [] = Putintofiles(obj, thisone, thisplot)
        type = '.png';
        tog = strcat(obj.newname, thisplot, type);
        slash = '/';
        loc = strcat(obj.imagedir, 'Images', slash, tog);
        saveas(thisone, loc);
        pause(obj.figs)
        close(thisone);
      end
   end
end
