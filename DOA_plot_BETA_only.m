clear all
close all
clc

addpath(genpath('../SDMtools/'));
addpath(genpath('../data/'));
addpath(genpath('functions/'));

%% load recorded data

% lade Daten mit mat file
filepath='../data/RIRs/';
%filepath='D:\Daten\DORAV\Parameter\SDM_M\';


Files=dir([filepath,'*.mat']);

if(length(Files)>1)
    for idx=1 : length(Files)
        disp([num2str(idx),': ',Files(idx).name])
    end
    prompt = 'select file number which should be loaded:';
    kselect=input(prompt);
else
    kselect=1;
end


filename=Files(kselect).name;
load([filepath filename])


for idxSpeaker=1 : 1 : length(irs.speakerNames)
    
    fs              = double(irs.fs);    % Sampling Rate (in Hz). Only 48 kHz is recommended. Other sampling rates have not been tested.
    
    %RAUMABHÄNGIG----------------------------------------------------------
    MixingTime      = 0.1;                 % Mixing time (in seconds) of the room for rendering. Data after the mixing time will be rendered
    BRIRLength      = 0.25;                  % Duration of the rendered BRIRs (in seconds)

    SpeedSound      = 345;                  % Speed of sound in m/s (for SDM Toolbox DOA analysis)
    
    ArrayGeometry = [1 0 0;
                     0 -0.7071 -0.7071;
                     0 -0.7071 0.7071;
                     -1 0 0;
                     0 0.7071 0.7071;
                     0 0.7071 -0.7071;
                     0 0 0] * (0.1/2);
    
    SDM_Struct = createSDMStruct('c',SpeedSound,...
        'fs',irs.fs,...ab
        'micLocs',ArrayGeometry,...
        'winLen',62);
     
    SRIR = irs.ir{idxSpeaker}(:,:);
    
    %f_gain = [5.899974163771157,5.960146874493056,3.669551814562449,4.178224276367109,3.710109648676412,5.230395430311390,1]; %Messung RMS von IR
    %SRIR = SRIR .* f_gain;
    


    DOA = SDMPar(SRIR, SDM_Struct); %Direction of arrival
    P = SRIR(:,7); %Impulse response
    


    %% from here on only plot building
    
    res = 1; % DOA resolution of the polar response
    t_steps = 1; %in ms
    
     
    ir_threshold = 0.5; % treshold level for the beginning of direct sound
     
    %Find the direct sound
    ind = find(abs(P)/max(abs(P)) > ir_threshold,1,'first');
    pre_threshold = round(0.001*fs); % go back 1 ms
    t_ds = ind - pre_threshold; % direct sound begins from this sample
    
    % make sure that the time index of direct sound is greater than or equal to 1
    t_ds = max(t_ds, 1);
    
    t_end = MixingTime/(t_steps/1000);
    
    ts = round( 1 : (t_steps/1000*fs) : (t_end/1000*fs)+(t_steps/1000*fs) );
    
    % Iterate through different time windows
    for k = 1:length(ts)-1
        
        t1 = t_ds+ts(k);
        t2 = t_ds+ts(k+1);
        tmpDOA = DOA(t1:t2,:);
        tmpP = P(t1:t2);
        
        %t2 = min(t_ds+t(k),(t_end/1000*fs));
        %tmpDOA = DOA{idx_s}(t2:(t_end/1000*fs),:);
        %tmpP = P{idx_s}(t2:(t_end/1000*fs));
        
        az_corr = 0;
        el_corr = 0;
        
        [az,el,~] = cart2sph(tmpDOA(:,1),tmpDOA(:,2),tmpDOA(:,3));
        [x,y,z] = sph2cart(az+az_corr,el+el_corr,1);
        [az,el,~] = cart2sph(x,y,z);
        
        % Find the closest direction in the grid for each image-source
        AZ = round(az(:)*180/pi/res)*res;
        EL = round(el(:)*180/pi/res)*res;
        
        % Pressure to energy
        A2 = tmpP.^2;
        A2 = A2(:);
        
        % Doughnut weighting for angles, i.e. cosine weighting
        doa_az_vec = -180:res:180;
        doa_el_vec = -90:res:90;
        
        H = zeros(length(doa_az_vec), length(doa_el_vec));
        
        for idx = 1 : length(A2)
            if( ~((isnan(AZ(idx))) && (isnan(EL(idx)))) )
                idxAz=find(AZ(idx)==doa_az_vec);
                idxEl=find(EL(idx)==doa_el_vec);
                H(idxAz,idxEl) = H(idxAz,idxEl) + A2(idx);
            end
        end
        HS(:,:,k) = H;
       
    end
    
    noisefloor_dB = -100;
    maxdB = 10*log10(max(max(max(HS))));
    %ges_dB = 10*log10(sum(HS,'all'));
    %mph = 0.0001*(10^(maxdB/10));
    mph = 10^(noisefloor_dB/10);
    
    
    az_tic =-180 :30: 180;
    el_tic =-90 :30: 90;
    
    Scale = 200;
    
    %Vorbereitung Plot DoA IR^2
    for idxk = k : -1 : 1
        
        tmpH =  HS(:,:,idxk);
                
        count=0;
        for idx = 1:size(tmpH,1)
            if(find(tmpH(idx,:)>mph,1))
                [pks,loc] = findpeaks(tmpH(idx,:),'MinPeakHeight',mph);
                %[pks,loc] = findpeaks(tmpH(idx,:));
                if(~isempty(pks))
                    for idxC = 1 : length(pks)
                        count=count+1;
                        %location az el
                        idxPeaksAzEl(count,:)=[idx,loc(idxC)];
                    end
                end
            end
        end
        
        for idx = 1 : count
            tmpA2(idx) = tmpH(idxPeaksAzEl(idx,1),idxPeaksAzEl(idx,2));
        end
        
        if(count~=0)
            peak_Az_El_R(idxk)={[idxPeaksAzEl,tmpA2']};  
        end
        clear A tmpA2 idxPeaksAzEl tmpH;
    end
    
    
    
    
    
    %light and dark mode
    
    plot_mode = "dark";
    %Sense_of_Rotation_Az = "clockwise"
    Sense_of_Rotation_Az = "counterclockwise";
    
    switch(plot_mode)
        case "dark"
            cmap = colormap('hot');
            stepw_color = floor(230/(length(ts)-1));
            %newcolors = cmap(1:stepw_color:256,:);
            newcolors = cmap(230:-stepw_color:1,:);
            BackColor = [0.1 0.1 0.1];
            GridColor = [0.9 0.9 0.9];
            AlphaColor = 0.8;
        case "light"
            cmap = colormap('jet');
            stepw_color = floor(230/(length(ts)-1));
            newcolors = cmap(1:stepw_color:230,:);
            %newcolors = cmap(230:-stepw_color:1,:);
            BackColor = [1 1 1];
            GridColor = [0.1 0.1 0.1];
            AlphaColor = 0.8;
        otherwise
    end
    
    if(length(newcolors)<length(ts))
       error('Farbabstufung zu klein') 
    end
    
    DefFonSiz = 14;
    %% Plot abs(IR)(abs pressure) over time
    figP=figure;
    set(figP,'defaultfigurecolor',[1 1 1],'DefaultAxesFontSize',DefFonSiz);
    for k = 1:length(ts)-1
        t1 = t_ds+ts(k);
        t2 = t_ds+ts(k+1);
        plot((ts(k):ts(k+1))/fs,(abs(P(t1:t2)/max(abs(P)))),'Color',[newcolors(k,:) AlphaColor])
        %plot((ts(k):ts(k+1))/fs,10*log10(abs(P{idx_s}(t1:t2)/max(abs(P{idx_s}))).^2),'Color',[newcolors(k,:) AlphaColor])
        %t2 = min(t_ds+t(k),(t_end/1000*fs));
        %plot((t2:(t_end/1000*fs)),P{idx_s}(t2:(t_end/1000*fs)),'Color',newcolors(k,:))
        hold on
    end
    grid on
    hold off
    xlim([0 t_end/1000])
    xlabel('time in s')
    ylabel('absolute amplitude')
    title(irs.room)
    
    ax = gca;
    ax.Color = BackColor;
    ax.GridColor = GridColor;
    close;
    
    
    %% Plot DoA IR (abs pressure)
    %DefFonSiz = 10;
        
    figP = figure;
    %set(figP, 'Position', [300, 150, 1024, 768],'defaultfigurecolor',[1 1 1],'DefaultAxesFontSize',DefFonSiz);
    %axesm ('apianus', 'Frame', 'on', 'Grid', 'on');
    set(figP, 'Position', [300, 150, 1024, 496],'defaultfigurecolor',[1 1 1],'DefaultAxesFontSize',DefFonSiz);
    
    %pic=imread('D:\Daten\PWHoertest_2021\H1562_E.jpg');
    %imagesc(180:-res:-180,-90:res:90,pic)
    %hold on
    
    %ax = axesm ('apianus', 'Frame', 'on', 'Grid', 'on');
    
    for idxk = length(peak_Az_El_R) : -1 : 1
        
        tmp = peak_Az_El_R{idxk};
        if(~isempty(tmp))
            
            A =sqrt(tmp(:,3)./((10^(maxdB/10))))*Scale;
            tmp=tmp(find(A~=0),:);
            
            scatter(doa_az_vec(tmp(:,1)),doa_el_vec(tmp(:,2)),A,'filled','MarkerFaceColor',newcolors(idxk,:),'MarkerFaceAlpha',AlphaColor)
            
            %scatterm(doa_az_vec(tmp(:,1)),doa_el_vec(tmp(:,2)),A,A,'filled','MarkerFaceColor',newcolors(idxk,:),'MarkerFaceAlpha',AlphaColor)
            hold on
        end
    end
    grid on
    xticks(az_tic)
    yticks(el_tic)
    xlim([-180 180])
    ylim([-90 90])
    title(irs.room)
    xlabel("azimuth in degree")
    ylabel("elevation in degree")
    switch(Sense_of_Rotation_Az)
        case 'clockwise'
            view([0 90])
        case 'counterclockwise'
            view([180 -90])
        otherwise
    end
    ax = gca;
    ax.Color = BackColor;
    ax.GridColor = GridColor;
    close;
    

    %% plot hist az und el over time

    tmp = sqrt(HS);

    [s_az, s_el, s_t] = size(HS);

    tmpAz = squeeze(sum(tmp,2));
    %tmpAz = tmpAz/max(max(tmpAz)); %Normiert auf direktschall
    tmpAz = tmpAz/max(sum(tmpAz,2)); %Normiert auf über die Zeit summierte absolute amplitude

    tmpEl = squeeze(sum(tmp,1));
    %tmpEl = tmpEl/max(max(tmpEl)); %Normiert auf direktschall
    tmpEl = tmpEl/max(sum(tmpEl,2)); %Normiert auf über die Zeit summierte absolute amplitude

    if(0)
        figP = figure;
        set(figP, 'Position', [300, 150, 1024, 496],'defaultfigurecolor',[1 1 1],'DefaultAxesFontSize',DefFonSiz);
        bp=bar(tmpAz,'stacked','FaceColor','flat');
        for k = 1:s_t
            bp(k).FaceColor = newcolors(k,:);
        end
        ylim([0 1])
        xticks(1 : (s_az-1)/(length(az_tic)-1) : s_az)
        xticklabels(az_tic)
        xlim([1 s_az])
        xlabel("azimuth in degree")
        ylabel("summed normalized absolute amplitude")
        switch(Sense_of_Rotation_Az)
            case 'clockwise'
                view([0 90])
            case 'counterclockwise'
                view([180 -90])
            otherwise
        end
        grid on
        ax = gca;
        ax.Color = BackColor;
        ax.GridColor = GridColor;
        close;

        figP = figure;
        set(figP, 'Position', [300, 150, 512, 496],'defaultfigurecolor',[1 1 1],'DefaultAxesFontSize',DefFonSiz);
        bp=bar(tmpEl,'stacked','FaceColor','flat');
        for k = 1:s_t
            bp(k).FaceColor = newcolors(k,:);
        end
        ylim([0 1])
        xticks(1 : (s_el-1)/(length(el_tic)-1) : s_el)
        xlim([1 s_el])
        xticklabels(el_tic)
        xlabel("elevation in degree")
        ylabel("summed normalized absolute amplitude")
        view([90 -90])
        grid on
        ax = gca;
        ax.Color = BackColor;
        ax.GridColor = GridColor;
        close;
    end
    
    %% plot All in one
    
    DefFonSiz=14;
    figP = figure;
    set(figP, 'Position', [100, 100, 1000, 500],'defaultfigurecolor',[1 1 1]);
    %set(figP, 'Position', [0, 0, 1500, 1100],'defaultfigurecolor',[1 1 1],'DefaultAxesFontSize',DefFonSiz);

    tiledplot=tiledlayout(2,3);
    tiledplot.TileSpacing = 'compact';
    %title(tiledplot,irs.room,'FontWeight','bold')
    title(tiledplot,['room: ',char(irs.room),'; speaker name: ',char(irs.speakerNames(idxSpeaker))],'FontWeight','bold','FontSize',DefFonSiz+2);
    
    nexttile
    % p over el
    bp=bar(tmpEl,'stacked','FaceColor','flat');
    for k = 1:s_t
        bp(k).FaceColor = newcolors(k,:);
    end
    ylim([0 1])
    xticks(1 : (s_el-1)/(length(el_tic)-1) : s_el)
    xlim([1 s_el])
    xticklabels(el_tic)
    xlabel("elevation in degree")
    ylabel("summed normalized absolute amplitude")
    view([90 -90])
    grid on
    ax = gca;
    ax.Color = BackColor;
    ax.GridColor = GridColor;
    
    nexttile([1 2])
    %Kreisplot az el
    for idxk = length(peak_Az_El_R) : -1 : 1
        
        tmp = peak_Az_El_R{idxk};
        if(~isempty(tmp))
            
            A =sqrt(tmp(:,3)./((10^(maxdB/10))))*Scale;
            tmp=tmp(find(A~=0),:);
            
            scatter(doa_az_vec(tmp(:,1)),doa_el_vec(tmp(:,2)),A,'filled','MarkerFaceColor',newcolors(idxk,:),'MarkerFaceAlpha',AlphaColor)
            hold on
        end
    end
    grid on
    xticks(az_tic)
    yticks(el_tic)
    xlim([-180 180])
    ylim([-90 90])
    xlabel("azimuth in degree")
    ylabel("elevation in degree")
    switch(Sense_of_Rotation_Az)
        case 'clockwise'
            view([0 90])
        case 'counterclockwise'
            view([180 -90])
        otherwise
    end
    ax = gca;
    ax.Color = BackColor;
    ax.GridColor = GridColor;
    
    nexttile
    %plot IR
    for k = 1:length(ts)-1
        t1 = t_ds+ts(k);
        t2 = t_ds+ts(k+1);
        plot((ts(k):ts(k+1))/fs,(abs(P(t1:t2)/max(abs(P)))),'Color',[newcolors(k,:) AlphaColor])
        %plot((ts(k):ts(k+1))/fs,10*log10(abs(P{idx_s}(t1:t2)/max(abs(P{idx_s}))).^2),'Color',[newcolors(k,:) AlphaColor])
        %t2 = min(t_ds+t(k),(t_end/1000*fs));
        %plot((t2:(t_end/1000*fs)),P{idx_s}(t2:(t_end/1000*fs)),'Color',newcolors(k,:))
        hold on
    end
    grid on
    hold off
    xlim([0 t_end/1000])
    xlabel('time in s')
    ylabel('absolute amplitude')
    ax = gca;
    ax.Color = BackColor;
    ax.GridColor = GridColor;
    
    nexttile([1 2])
    
    bp=bar(tmpAz,'stacked','FaceColor','flat');
    for k = 1:s_t
        bp(k).FaceColor = newcolors(k,:);
    end
    ylim([0 1])
    xticks(1 : (s_az-1)/(length(az_tic)-1) : s_az)
    xticklabels(az_tic)
    xlim([1 s_az])
    xlabel("azimuth in degree")
    ylabel("summed normalized absolute amplitude")
    switch(Sense_of_Rotation_Az)
        case 'clockwise'
            view([0 90])
        case 'counterclockwise'
            view([180 -90])
        otherwise
    end
    grid on
    ax = gca;
    ax.Color = BackColor;
    ax.GridColor = GridColor;
    
    
    
end

