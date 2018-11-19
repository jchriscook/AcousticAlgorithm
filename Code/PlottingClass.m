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
      % the ranged vals
      Rset
   end
   
   methods
       function setzone = setting(obj, data)
            % create a time array
           time = linspace(0, 20, 1200000);
           
           if obj.Range(1) == 0
               begin = 1;
           else 
               begin = obj.Range(1);
           end
           
           if obj.Range(2) ~= 1
               close = obj.Range(2);
               stop = close * length(time);
           end
           start = begin * length(time);
  
           if exist('stop') == 1
               zz = data(start:stop);
           else
               zz = data(start:end);
           end
           
           setzone = zz;
       end
       function [] = Perform(obj)
           %% Fourier Domain - Cook and Weber      
           L = length(obj.R1);
           % the length will always have a certain delta t associated
           timestep = 1200 / L;
           % the start and the stop view of the data given the range 
           % ----- Finish Later ----- %
           % beginning = obj.range(1) * L;
           % stop = obj.range(2) * L;
           % ----- 
           
           % sample rate 
           fs = 1 / timestep;        % Hz !
           FoldingFreq = fs / 2;
           
           thisone = figure;
           set(thisone, 'Visible', obj.figs);
           
          
           
           %% Linear Magnitude - Subsection (C&W)
           subplot(2,2,1)
           % pass the object into setting method
           variable = obj.setting(obj.R1);
           
           Yx = fft(variable);
           n1 = length(variable)/2;
           n1 = ceil(n1);
           freq1 = linspace(0, FoldingFreq, n1);
           amp_specx = abs(Yx)/n1;
           db_spec1 = mag2db(amp_specx);
           stem(freq1,amp_specx(1:n1), '-b');
           xlim([0 100]);
           ylim([0 0.05]);
           title('Microphone 1: Roof','fontsize',15);
           xlabel('Frequency (Hz)','fontsize',15);
           ylabel('Linear Magnitude','fontsize',15);
           %set(gca,'Ydir','reverse')
           set(gca,'FontSize',12);
      
           subplot(2,2,2)
           % pass the object into setting method
           variable2 = obj.setting(obj.R2);
           
           Yy = fft(variable2);
           n2 = length(variable2)/2;
           n2 = ceil(n2);
           freq2 = linspace(0,FoldingFreq,n2);
           amp_specy = abs(Yy)/n2;
           db_spec2 = mag2db(amp_specy);
           stem(freq2,amp_specy(1:n2), '-g');
           xlim([0 100]);
           ylim([0 0.05]);
           title('Microphone 2: South','fontsize',15);
           xlabel('Frequency (Hz)','fontsize',15);
           ylabel('Linear Magnitude','fontsize',15);
           %set(gca,'Ydir','reverse')
           set(gca,'FontSize',12);

           subplot(2,2,[3, 4])
           % pass the object into setting method
           variable3 = obj.setting(obj.R3);
           
           Yz = fft(variable3);
           n3 = length(variable3)/2;
           n3 = ceil(n3);
           freq3 = linspace(0,FoldingFreq,n3);
           amp_specz = abs(Yz)/n3;
           db_spec3 = mag2db(amp_specz);
           stem(freq3,amp_specz(1:n3), '-r');
           xlim([0 100]);
           ylim([0 0.05]);
           title('Microphone 3: North','fontsize',15);
           xlabel('Frequency (Hz)','fontsize',15);
           ylabel('Linear Magnitude','fontsize',15);
           %set(gca,'Ydir','reverse')
           set(gca,'FontSize',12);
           
           % save the figure
           thisplot = ' - freqdomain';
           obj.Putintofiles(thisone, thisplot);          
           %% Decibal Magnitude
           thisone2 = figure;
           set(thisone2, 'Visible', obj.figs);
           
           
           subplot(2,2,1)
           plot(freq1, db_spec1(1:n1), '-b');
           ylim([-200 0]);
           xlim([0 100]);
           title('Microphone 1: Roof','fontsize',15);
           xlabel('Frequency (Hz)','fontsize',15);
           ylabel('Decibal Magnitude (dB)','fontsize',15);
           set(gca,'FontSize',12);
      
           subplot(2,2,2)
           plot(freq2, db_spec2(1:n2), '-g');
           ylim([-200 0]);
           xlim([0 100]);
           title('Microphone 2: South','fontsize',15);
           xlabel('Frequency (Hz)','fontsize',15);
           ylabel('Decibal Magnitude (dB)','fontsize',15);
           set(gca,'FontSize',12);

           subplot(2,2,[3, 4])
           plot(freq3, db_spec3(1:n3), '-r');
           ylim([-200 0]);
           xlim([0 100]);
           title('Microphone 3: North','fontsize',15);
           xlabel('Frequency (Hz)','fontsize',15);
           ylabel('Decibal Magnitude (dB)','fontsize',15);
           set(gca,'FontSize',12);
           
           thisplot = ' - db';
           obj.Putintofiles(thisone2, thisplot);         
       end
       
      function [] = filtering(obj, s, name)
          %% Lowpass filter (Butterworth) - Elbing
          
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
