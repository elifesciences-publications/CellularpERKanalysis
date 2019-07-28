%% AGGREGATE ALL cellular pERK data. 
% This code is specifically meant for the analysis of files within the downloadble folder "processed data for upload"
% This code will require modification for other types of data 
% resampleindices.m is a dependency 

clf; clear all;

close all;
%change folder path accordingly - where you want the output to be saved
savefolder = '/Users/carolinewee/Desktop/pERK/epsfiles/'

%change folder path accordingly - path to the the data folder (each
%stimulus should belong to a different subfolder)
folderpath = '/Users/carolinewee/Desktop/pERK/Quantitative/processeddataforupload/'
folders = dir(folderpath)

% generate a structure array to store all the data

field1 = 'Shock_60V', value1 = [];
field2 = 'MO_10uM', value2 = [];
field3 = 'Heat_37deg', value3 = [];
field4 = 'Taps', value4 = [];
field5 = 'MO_25uM', value5 = [];
field6 = 'Shock_50V', value6 = [];
field7 = 'Shock_40V', value7 = [];
field8 = 'Heat_40deg', value8 = [];

alldata = struct(field1, value1, field2, value2, field3, value3, field4, value4, field5, value5, field6, value6, field7, value7, field8, value8);

subfield1 = 'data_control';
subfield2 = 'data_control_FLUOR';
subfield3 = 'data_treatment';
subfield4 = 'data_treatment_FLUOR';

fnames = fieldnames(alldata);
for i = 1:length(fieldnames(alldata))
    
    exp_type = char(fnames(i));
    alldata.(exp_type) = struct(subfield1, [], subfield2, [], subfield3, [], subfield4, []);
    
end

%now loop through all folders (MO, Shock, Taps, Heat etc)
for i= 1: length(folders)
    
    name = folders(i).name
    
    findstring = strfind(name, '.');
    if isempty(findstring)
    else
        continue
    end
    
    files = dir(strcat(folderpath, '/', name));
    
    %loop through all csv files in the folder that contains the data
    for j =1:length(files)
        
        filename = files(j).name;
        findstring2 = strfind(filename, '.csv');
        
        if ~isempty(findstring2)
            rawdata = load(sprintf('%s%s/%s', folderpath, name, filename));
        else continue
        end
        
        findstring_control = strfind(filename, 'control');
        if ~isempty(findstring_control)
            condition = 'control';
        else
            condition = 'treatment';
        end
        
        %now to identify if GFP or other channel
        findstring_GFP = strfind(filename, 'GFP');
        findstring_allcells = strfind(filename, 'allcells');
        
        if ~isempty(findstring_GFP)
            celltype = 'GFP';
        elseif ~isempty(findstring_allcells)
            celltype = 'allcells';
        end
        
        datatype = strcat(condition, celltype);
        
        %my naming convention
        date = filename(1:8);
        datename = strcat('exp', date);
        
        if strcmp(condition, 'control') && strcmp(celltype, 'allcells')
            subfield = 'data_control';
        elseif strcmp(condition, 'control') && strcmp(celltype, 'GFP')
            subfield = 'data_control_FLUOR';
        elseif strcmp(condition, 'treatment') && strcmp(celltype, 'allcells')
            subfield = 'data_treatment';
        elseif strcmp(condition, 'treatment') && strcmp(celltype, 'GFP')
            subfield = 'data_treatment_FLUOR';
        end
        
        %this is either to set field (i.e. date), or add a new field if
        %it is a new date.
        alldata.(name).(subfield) = setfield(alldata.(name).(subfield), datename, rawdata);
        
    end
end

cd('/Users/carolinewee/Desktop/pERK/Quantitative/datafiles/')
save('alldissectedpERKdata.mat', 'alldata')

%% ANALYZE AND AGGREGATE DATA
clear rawdata
% if running from this section, may need to restate savefolder below
%change folder path accordingly
%savefolder = '/Users/carolinewee/Desktop/pERK/epsfiles/'
cd(savefolder);

% Currently for 8 fields
analyzeddata = struct(field1, value1, field2, value2, field3, value3, field4, value4, field5, value5, field6, value6, field7, value7, field8, value8);
fnames_exp = fieldnames(alldata);

for e = 1:length(fnames_exp)
    
    exp_type = char(fnames_exp(e));
    
    % modifiable -- placed here so won't be cleared
    useERK = 0;
    f = 18;
    oxtcontrol_color = 'k'
    oxttreatment_color = [204/255 0 0];
    otherscontrol_color = [0.5 0.5 0.5];
    otherstreatment_color = [0 0 102/255];
    
    % to concatenate fish from all the dates
    fnames = fieldnames(alldata.(exp_type).data_control);
    
    for d = 1:length(fnames)
        
        exp_date = char(fnames(d));  % to specify which date to analyze
        
        % to set up analyzed data array
        analyzeddata.(exp_type) = setfield(analyzeddata.(exp_type), exp_date, []);
        
        % to load data
        data_treatment = alldata.(exp_type).data_treatment.(exp_date);
        data_treatment_FLUOR = alldata.(exp_type).data_treatment_FLUOR.(char(exp_date));
        data_control = alldata.(exp_type).data_control.(exp_date);
        data_control_FLUOR= alldata.(exp_type).data_control_FLUOR.(exp_date);
        
        % analyze treatment data
        ratios_treatment = [];
        GFP_treatment = [];
        Xpos_treatment = [];
        Ypos_treatment = [];
        Zpos_treatment = [];
        
        FISH_fluorescence_treatment = data_treatment(:, 1:5:end);
        FISH_xpos_treatment = data_treatment(:, 2:5:end);
        FISH_ypos_treatment = data_treatment(:, 3:5:end);
        FISH_zpos_treatment = data_treatment(:, 5:5:end);
        
        for i = 1: size(FISH_fluorescence_treatment,2)
            
            for j = 2:3:size(FISH_fluorescence_treatment,1)
                
                if useERK == 1
                    ratios_treatment((j-2)/3+1,i) = FISH_fluorescence_treatment(j,i)/FISH_fluorescence_treatment(j+1,i);
                else
                    ratios_treatment((j-2)/3+1,i) = FISH_fluorescence_treatment(j,i);
                end
                
                GFP_treatment((j-2)/3+1,i) = FISH_fluorescence_treatment(j-1,i);
                Xpos_treatment((j-2)/3+1,i) = FISH_xpos_treatment(j-1,i);
                Ypos_treatment((j-2)/3+1,i) = FISH_ypos_treatment(j-1,i);
                Zpos_treatment((j-2)/3+1,i) = FISH_zpos_treatment(j-1,i);
            end
        end
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_treatment', ratios_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'GFP_treatment', GFP_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Xpos_treatment', Xpos_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Ypos_treatment', Ypos_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Zpos_treatment', Zpos_treatment);
        
        % analyze control data
        ratios_control = [];
        GFP_control  = [];
        Xpos_control  = [];
        Ypos_control  = [];
        Zpos_control  = [];
        
        FISH_fluorescence_control  = data_control(:, 1:5:end);
        FISH_xpos_control  = data_control(:, 2:5:end);
        FISH_ypos_control  = data_control(:, 3:5:end);
        FISH_zpos_control  = data_control(:, 5:5:end);
        
        
        for i = 1:size(FISH_fluorescence_control ,2)
            
            for j = 2:3:size(FISH_fluorescence_control,1)
                if useERK == 1
                    ratios_control((j-2)/3+1,i) = FISH_fluorescence_control (j,i)/FISH_fluorescence_control(j+1,i);
                else
                    ratios_control((j-2)/3+1,i) = FISH_fluorescence_control (j,i);
                end
                
                GFP_control((j-2)/3+1,i) = FISH_fluorescence_control (j-1,i);
                Xpos_control((j-2)/3+1,i) = FISH_xpos_control (j-1,i);
                Ypos_control((j-2)/3+1,i) = FISH_ypos_control (j-1,i);
                Zpos_control((j-2)/3+1,i) = FISH_zpos_control (j-1,i);
            end
        end
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_control', ratios_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'GFP_control', GFP_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Xpos_control', Xpos_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Ypos_control', Ypos_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Zpos_control', Zpos_control);
        
        % Consolidate pERK/tERK ratios and GFP intensities into single array
        controldata = ratios_control(~isnan(ratios_control));
        treatmentdata = ratios_treatment(~isnan(ratios_treatment));
        controlGFP = GFP_control(~isnan(GFP_control));
        treatmentGFP = GFP_treatment(~isnan(GFP_treatment));
        
        % this is to consolidate xyzpositions into single array
        normYPOS_control = [];
        normYPOS_treatment = [];
        normZPOS_control = [];
        normZPOS_treatment = [];
        normXPOS_control = [];
        normXPOS_treatment = [];
        
        for i = 1:size(Ypos_control,2)
            normYPOS_control(:,i)=Ypos_control(:,i)-min(Ypos_control(:,i)); %normalize to most anterior
            normZPOS_control(:,i)=Zpos_control(:,i)-min(Zpos_control(:,i)); %normalize to most dorsal (or ventral if mounted upside down).
            normXPOS_control(:,i)=Xpos_control(:,i)-min(Xpos_control(:,i));
        end
        
        for i = 1:size(Ypos_treatment,2)
            normYPOS_treatment(:,i)=Ypos_treatment(:,i)-min(Ypos_treatment(:,i)); %normalize to most anterior
            normZPOS_treatment(:,i)=Zpos_treatment(:,i)-min(Zpos_treatment(:,i));
            normXPOS_treatment(:,i)=Xpos_treatment(:,i)-min(Xpos_treatment(:,i));
        end
        
        yposdata_control = normYPOS_control(~isnan(normYPOS_control));
        yposdata_treatment= normYPOS_treatment(~isnan(normYPOS_treatment));
        zposdata_control = normZPOS_control(~isnan(normZPOS_control));
        zposdata_treatment= normZPOS_treatment(~isnan(normZPOS_treatment));
        xposdata_control = normXPOS_control(~isnan(normXPOS_control));
        xposdata_treatment= normXPOS_treatment(~isnan(normXPOS_treatment));
        
        % Extract OXT (GFP intensity >criterion) neurons from tERK image
        oxtcriterion = 400; %% this can be changed
        idx_oxt_control = find(controlGFP > oxtcriterion);
        idx_oxt_treatment = find(treatmentGFP > oxtcriterion);
        
        normfactor = nanmean(controldata);
        controldata = controldata/normfactor;  %normalize control data
        all_controldata = controldata; %store additional copy as another variable
        treatmentdata = treatmentdata/normfactor; %normalize treatment data
        all_treatmentdata = treatmentdata;
        
        oxt_controldata = all_controldata(idx_oxt_control);
        oxt_treatmentdata = all_treatmentdata(idx_oxt_treatment);
        
        % extract GFP-positive cells from ratios_treatment and ratios_control array to preserve fish information
        ratios_control_others = (GFP_control< oxtcriterion).*ratios_control;
        ratios_control_others(ratios_control_others ==0) = nan;
        ratios_control_oxt = (GFP_control>oxtcriterion).*ratios_control;
        ratios_control_oxt(ratios_control_oxt==0) = nan;
        
        ratios_treatment_others = (GFP_treatment<oxtcriterion).*ratios_treatment;
        ratios_treatment_others(ratios_treatment_others ==0) = nan;
        ratios_treatment_oxt = (GFP_treatment>oxtcriterion).*ratios_treatment;
        ratios_treatment_oxt(ratios_treatment_oxt ==0) = nan;
        
        % now normalize before saving
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_control_oxt', ratios_control_oxt/normfactor);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_control_others', ratios_control_others/normfactor);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_treatment_oxt', ratios_treatment_oxt/normfactor);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_treatment_others', ratios_treatment_others/normfactor);
        
        figure(1);
        clf;
        [f_cont,xi_cont] = ksdensity(oxt_controldata);
        [f_treatment,xi] = ksdensity(oxt_treatmentdata);
        
        plot(xi_cont, f_cont, 'Color', 'k', 'LineWidth', 2);
        hold on;
        plot(xi, f_treatment, 'Color', 'r', 'LineWidth', 2);
        set(gca, 'FontSize', f)
        box on;
        set(gca, 'FontSize', f)
        set(gca,'YTick',[])
        set(gca,'YTicklabel','')
        box off;
        hold on;
        
        [h_ks, p_ks] = kstest2(oxt_controldata, oxt_treatmentdata)
        [p_rs, h_rs] = ranksum(oxt_controldata, oxt_treatmentdata)
        
        % Now get rid of oxt data to analyze data from other neurons
        controldata(idx_oxt_control) = nan;
        treatmentdata(idx_oxt_treatment) = nan;
        others_controldata = controldata(~isnan(controldata));
        others_treatmentdata = treatmentdata(~isnan(treatmentdata));
        
        figure(1);
        [f_cont,xi_cont] = ksdensity(others_controldata);
        [f_treatment,xi] = ksdensity(others_treatmentdata);
        
        plot(xi_cont, f_cont, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
        hold on;
        plot(xi, f_treatment, 'Color', 'b', 'LineWidth', 2);
        set(gca, 'FontSize', f)
        box on;
        set(gca, 'FontSize', f)
        set(gca,'YTick',[])
        set(gca,'YTicklabel','')
        %set(gca,'XTick',[])
        %set(gca,'XTicklabel','')
        axis([0 2.5 0 4])
        box off;
        
        print(sprintf('%s_%s%s', exp_type, exp_date,'-allcellshistogram'), '-depsc2');
        
        [h_ks2, p_ks2] = kstest2(others_controldata, others_treatmentdata)
        [p_rs2, h_rs2] = ranksum(others_controldata, others_treatmentdata)
        
        % now save data extracted from all cells
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_controldata', oxt_controldata);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_treatmentdata', oxt_treatmentdata);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_controldata', others_controldata);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_treatmentdata', others_treatmentdata);
        
        % do the same for x y and zpos
        oxt_yposdata_control =yposdata_control(idx_oxt_control);
        oxt_yposdata_treatment = yposdata_treatment(idx_oxt_treatment);
        oxt_zposdata_control =zposdata_control(idx_oxt_control);
        oxt_zposdata_treatment = zposdata_treatment(idx_oxt_treatment);
        oxt_xposdata_control =xposdata_control(idx_oxt_control);
        oxt_xposdata_treatment = xposdata_treatment(idx_oxt_treatment);
        
        % new vector to store all data
        all_yposdata_control = yposdata_control;
        all_yposdata_treatment = yposdata_treatment;
        all_zposdata_control = zposdata_control;
        all_zposdata_treatment = zposdata_treatment;
        all_xposdata_control = xposdata_control;
        all_xposdata_treatment = xposdata_treatment;
        
        % now overwriting
        yposdata_control(idx_oxt_control) = nan;
        yposdata_treatment(idx_oxt_treatment) = nan;
        others_yposdata_control = yposdata_control(~isnan(yposdata_control));
        others_yposdata_treatment = yposdata_treatment(~isnan(yposdata_treatment));
        
        zposdata_control(idx_oxt_control) = nan;
        zposdata_treatment(idx_oxt_treatment) = nan;
        others_zposdata_control = zposdata_control(~isnan(zposdata_control));
        others_zposdata_treatment = zposdata_treatment(~isnan(zposdata_treatment));
        
        xposdata_control(idx_oxt_control) = nan;
        xposdata_treatment(idx_oxt_treatment) = nan;
        others_xposdata_control = xposdata_control(~isnan(xposdata_control));
        others_xposdata_treatment = xposdata_treatment(~isnan(xposdata_treatment));
        
        % now save data extracted from all cells
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_yposdata_control', oxt_yposdata_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_xposdata_control', oxt_xposdata_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_zposdata_control', oxt_zposdata_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_yposdata_control', others_yposdata_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_xposdata_control', others_xposdata_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_zposdata_control', others_zposdata_control);
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_yposdata_treatment', oxt_yposdata_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_xposdata_treatment', oxt_xposdata_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_zposdata_treatment', oxt_zposdata_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_yposdata_treatment', others_yposdata_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_xposdata_treatment', others_xposdata_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'others_zposdata_treatment', others_zposdata_treatment);
        
       
        %% NOW Use transgenic line data (where cells are extracted from fluorescence channel) for more accurate GFP reading
        
        % analyze treatment data
        ratios_treatment_FLUOR = [];
        GFP_treatment_FLUOR = [];
        Xpos_treatment_FLUOR = [];
        Ypos_treatment_FLUOR = [];
        Zpos_treatment_FLUOR = [];
        
        FISH_fluorescence_treatment_FLUOR = data_treatment_FLUOR(:, 1:5:end);
        FISH_xpos_treatment_FLUOR = data_treatment_FLUOR(:, 2:5:end);
        FISH_ypos_treatment_FLUOR = data_treatment_FLUOR(:, 3:5:end);
        FISH_zpos_treatment_FLUOR = data_treatment_FLUOR(:, 5:5:end);
        
        for i = 1:size(FISH_fluorescence_treatment_FLUOR,2)
            
            for j = 2:3:size(FISH_fluorescence_treatment_FLUOR,1)
                if useERK == 1
                    ratios_treatment_FLUOR((j-2)/3+1,i) = FISH_fluorescence_treatment_FLUOR(j,i)./FISH_fluorescence_treatment_FLUOR(j+1,i);
                else
                    ratios_treatment_FLUOR((j-2)/3+1,i) = FISH_fluorescence_treatment_FLUOR(j,i);
                end
                GFP_treatment_FLUOR((j-2)/3+1,i) = FISH_fluorescence_treatment_FLUOR(j-1,i);
                Xpos_treatment_FLUOR((j-2)/3+1,i) = FISH_xpos_treatment_FLUOR(j-1,i);
                Ypos_treatment_FLUOR((j-2)/3+1,i) = FISH_ypos_treatment_FLUOR(j-1,i);
                Zpos_treatment_FLUOR((j-2)/3+1,i) = FISH_zpos_treatment_FLUOR(j-1,i);
            end
        end
        
        % analyze control data
        ratios_control_FLUOR = [];
        GFP_control_FLUOR  = [];
        Xpos_control_FLUOR  = [];
        Ypos_control_FLUOR  = [];
        Zpos_control_FLUOR  = [];
        
        FISH_fluorescence_control_FLUOR = data_control_FLUOR(:, 1:5:end);
        FISH_xpos_control_FLUOR  = data_control_FLUOR(:, 2:5:end);
        FISH_ypos_control_FLUOR  = data_control_FLUOR(:, 3:5:end);
        FISH_zpos_control_FLUOR  = data_control_FLUOR(:, 5:5:end);
        
        for i = 1:size(FISH_fluorescence_control_FLUOR ,2)
            
            for j = 2:3:size(FISH_fluorescence_control_FLUOR,1)
                if useERK == 1
                    ratios_control_FLUOR((j-2)/3+1,i) = FISH_fluorescence_control_FLUOR(j,i)/FISH_fluorescence_control_FLUOR(j+1,i);
                else
                    ratios_control_FLUOR((j-2)/3+1,i) = FISH_fluorescence_control_FLUOR(j,i);
                end
                GFP_control_FLUOR ((j-2)/3+1,i) = FISH_fluorescence_control_FLUOR(j-1,i);
                Xpos_control_FLUOR((j-2)/3+1,i) = FISH_xpos_control_FLUOR(j-1,i);
                Ypos_control_FLUOR((j-2)/3+1,i) = FISH_ypos_control_FLUOR(j-1,i);
                Zpos_control_FLUOR((j-2)/3+1,i) = FISH_zpos_control_FLUOR(j-1,i);
            end
        end
        
        normfactor2 = nanmean(ratios_control_FLUOR(:)); %this is the same norm factor as later on, just that I want to normalize right now
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_treatment_FLUOR', ratios_treatment_FLUOR/normfactor2);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'GFP_treatment_FLUOR', GFP_treatment_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Xpos_treatment_FLUOR', Xpos_treatment_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Ypos_treatment_FLUOR', Ypos_treatment_FLUOR);
        analyzeddata.(exp_type).(exp_date)= setfield(analyzeddata.(exp_type).(exp_date), 'Zpos_treatment_FLUOR', Zpos_treatment_FLUOR);
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'ratios_control_FLUOR', ratios_control_FLUOR/normfactor2);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'GFP_control_FLUOR', GFP_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Xpos_control_FLUOR', Xpos_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Ypos_control_FLUOR', Ypos_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'Zpos_control_FLUOR', Zpos_control_FLUOR);
        
        % Consolidate pERK/tERK ratios and GFP intensities into single array
        controldata_FLUOR = ratios_control_FLUOR(~isnan(ratios_control_FLUOR));
        treatmentdata_FLUOR = ratios_treatment_FLUOR(~isnan(ratios_treatment_FLUOR));
        controlGFP_FLUOR = GFP_control_FLUOR(~isnan(GFP_control_FLUOR));
        treatmentGFP_FLUOR= GFP_treatment_FLUOR(~isnan(GFP_treatment_FLUOR));
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'controldata_FLUOR', controldata_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'treatmentdata_FLUOR', treatmentdata_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'controlGFP_FLUOR', controlGFP_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'treatmentGFP_FLUOR', treatmentGFP_FLUOR);
        
        % to consolidate xyzpositions into single array
        normYPOS_control_FLUOR = [];
        normYPOS_treatment_FLUOR = [];
        normZPOS_control_FLUOR = [];
        normZPOS_treatment_FLUOR = [];
        normXPOS_control_FLUOR = [];
        normXPOS_treatment_FLUOR = [];
        
        for i = 1:size(Ypos_control_FLUOR,2)
            normYPOS_control_FLUOR(:,i)=Ypos_control_FLUOR(:,i)-min(Ypos_control_FLUOR(:,i)); %normalize to most anterior
            normZPOS_control_FLUOR(:,i)=Zpos_control_FLUOR(:,i)-min(Zpos_control_FLUOR(:,i)); %normalize to most dorsal (or ventral if mounted upside down).
            normXPOS_control_FLUOR(:,i)=Xpos_control_FLUOR(:,i)-min(Xpos_control_FLUOR(:,i));
        end
        
        for i = 1:size(Ypos_treatment_FLUOR,2)
            normYPOS_treatment_FLUOR(:,i)=Ypos_treatment_FLUOR(:,i)-min(Ypos_treatment_FLUOR(:,i)); %normalize to most anterior
            normZPOS_treatment_FLUOR(:,i)=Zpos_treatment_FLUOR(:,i)-min(Zpos_treatment_FLUOR(:,i));
            normXPOS_treatment_FLUOR(:,i)=Xpos_treatment_FLUOR(:,i)-min(Xpos_treatment_FLUOR(:,i));
        end
        
        yposdata_control_FLUOR = normYPOS_control_FLUOR(~isnan(normYPOS_control_FLUOR));
        yposdata_treatment_FLUOR= normYPOS_treatment_FLUOR(~isnan(normYPOS_treatment_FLUOR));
        zposdata_control_FLUOR = normZPOS_control_FLUOR(~isnan(normZPOS_control_FLUOR));
        zposdata_treatment_FLUOR= normZPOS_treatment_FLUOR(~isnan(normZPOS_treatment_FLUOR));
        xposdata_control_FLUOR = normXPOS_control_FLUOR(~isnan(normXPOS_control_FLUOR));
        xposdata_treatment_FLUOR= normXPOS_treatment_FLUOR(~isnan(normXPOS_treatment_FLUOR));
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'yposdata_control_FLUOR', yposdata_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'yposdata_treatment_FLUOR', yposdata_treatment_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'zposdata_control_FLUOR', zposdata_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'zposdata_treatment_FLUOR', zposdata_treatment_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'xposdata_control_FLUOR', xposdata_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'xposdata_treatment_FLUOR', xposdata_treatment_FLUOR);
        
        % Normalize to control
        normfactor2 = nanmean(controldata_FLUOR);
        oxt_controldata_FLUOR = controldata_FLUOR./normfactor2;
        oxt_treatmentdata_FLUOR = treatmentdata_FLUOR./normfactor2;
        
        % save data
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_controldata_FLUOR', oxt_controldata_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'oxt_treatmentdata_FLUOR', oxt_treatmentdata_FLUOR);
        
        figure(2);
        clf;
        [f_cont,xi_cont] = ksdensity(oxt_controldata_FLUOR);
        [f_treatment,xi] = ksdensity(oxt_treatmentdata_FLUOR);
        
        %f_cont = pdf(contdata, xi_cont);
        %f_treatment = pdf(treatmentdata, xi);
        
        plot(xi_cont, f_cont, 'Color', 'k', 'LineWidth', 2);
        hold on;
        plot(xi, f_treatment, 'Color', 'r', 'LineWidth', 2);
        set(gca, 'FontSize', f)
        box on;
        set(gca, 'FontSize', f)
        set(gca,'YTick',[])
        set(gca,'YTicklabel','')
        %set(gca,'XTick',[])
        %set(gca,'XTicklabel','')
        box off;
        hold on;
        
        % this is the same data from previous section
        [f_cont,xi_cont] = ksdensity(others_controldata);
        [f_treatment,xi] = ksdensity(others_treatmentdata);
        
        %f_cont = pdf(contdata, xi_cont);
        %f_treatment = pdf(treatmentdata, xi);
        
        plot(xi_cont, f_cont, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
        hold on;
        plot(xi, f_treatment, 'Color', 'b', 'LineWidth', 2);
        set(gca, 'FontSize', f)
        box on;
        set(gca, 'FontSize', f)
        set(gca,'YTick',[])
        set(gca,'YTicklabel','')
        %set(gca,'XTick',[])
        %set(gca,'XTicklabel','')
        axis([0 2.5 0 4])
        box off;
        
        print(sprintf('%s_%s%s', exp_type, exp_date,'-allcellsandGFPhistogram'), '-depsc2');
        
        % percentage active and spatial distribution using GFP fluorescence

        % modify percentile here if necessary
        percentile = 3;
        active_threshold2 = findthreshold(oxt_controldata_FLUOR, percentile, 'above');
        
        idx_active_control_FLUOR = find(oxt_controldata_FLUOR>active_threshold2);
        idx_active_treatment_FLUOR = find(oxt_treatmentdata_FLUOR>active_threshold2);
        fraction_responsive_oxt_control_FLUOR = length(idx_active_control_FLUOR)/length(oxt_controldata_FLUOR);
        fraction_responsive_oxt_treatment_FLUOR = length(idx_active_treatment_FLUOR)/length(oxt_treatmentdata_FLUOR);
        
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'active_threshold2', active_threshold2);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'fraction_responsive_oxt_control_FLUOR', fraction_responsive_oxt_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'fraction_responsive_oxt_treatment_FLUOR', fraction_responsive_oxt_treatment_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'idx_active_control_FLUOR', idx_active_control_FLUOR);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'idx_active_treatment_FLUOR', idx_active_treatment_FLUOR);
        
        %% now use active threshold from FLUOR to calculate active neurons in other conditions
        
        % commented out as not relevant unless you want to stick to the 3% threshold rather than use the threshold defined in previous section
        % percentile = 3;
        % active_threshold = findthreshold(all_controldata, percentile, 'above'); 
        
        idx_active_control_others = find(others_controldata > active_threshold2);
        idx_active_treatment_others = find(others_treatmentdata > active_threshold2);
        idx_active_control_oxt = find(oxt_controldata > active_threshold2);
        idx_active_treatment_oxt = find(oxt_treatmentdata > active_threshold2);
        
        fraction_responsive_others_control =  length(idx_active_control_others)/length(others_controldata);
        fraction_responsive_oxt_control = length(idx_active_control_oxt)/length(oxt_controldata);
        fraction_responsive_others_treatment = length(idx_active_treatment_others)/length(others_treatmentdata);
        fraction_responsive_oxt_treatment = length(idx_active_treatment_oxt)/length(oxt_treatmentdata);
        
        %save data
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'active_threshold', active_threshold);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'fraction_responsive_others_control', fraction_responsive_others_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'fraction_responsive_oxt_control', fraction_responsive_oxt_control);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'fraction_responsive_others_treatment', fraction_responsive_others_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'fraction_responsive_oxt_treatment', fraction_responsive_oxt_treatment);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'idx_active_control_oxt', idx_active_control_oxt);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'idx_active_treatment_oxt',idx_active_treatment_oxt);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'idx_active_control_others', idx_active_control_others);
        analyzeddata.(exp_type).(exp_date) = setfield(analyzeddata.(exp_type).(exp_date), 'idx_active_treatment_others',idx_active_treatment_others);
        
        %plot spatial for each experiment and save
        figure(3);
        clf
        subplot(1,2,2); %treatment
        scatter(xposdata_treatment_FLUOR, yposdata_treatment_FLUOR,  50, [0.5 0.5 0.5], 'filled');
        hold on;
        scatter(xposdata_treatment_FLUOR(idx_active_treatment_FLUOR), yposdata_treatment_FLUOR(idx_active_treatment_FLUOR), 50, oxttreatment_color, 'filled');
        set(gca, 'Ydir', 'reverse');
        axis([-10 75 0 160])
        box off;
        axis off;
        
        subplot(1,2,1); %control
        scatter(xposdata_control_FLUOR, yposdata_control_FLUOR,  50, [0.5 0.5 0.5], 'filled');
        hold on;
        scatter(xposdata_control_FLUOR(idx_active_control_FLUOR), yposdata_control_FLUOR(idx_active_control_FLUOR), 50, oxttreatment_color, 'filled');
        set(gca, 'Ydir', 'reverse');
        axis([-10 75 0 160])
        box off;
        axis off;
        
        print(sprintf('%s_%s%s', exp_type, exp_date, '-spatial_oxt_fromGFP'), '-depsc2');
        
        % plot spatial oxt from all cells
        figure(4);
        clf;
        subplot(1,2,2); %treatment
        scatter(oxt_xposdata_treatment, oxt_yposdata_treatment,  50, [0.5 0.5 0.5], 'filled');
        hold on;
        scatter(oxt_xposdata_treatment(idx_active_treatment_oxt), oxt_yposdata_treatment(idx_active_treatment_oxt), 50, oxttreatment_color, 'filled');
        set(gca, 'Ydir', 'reverse');
        axis([-10 75 0 160])
        box off;
        axis off
        
        subplot(1,2,1);
        scatter(oxt_xposdata_control, oxt_yposdata_control,  50, [0.5 0.5 0.5], 'filled');
        hold on;
        scatter(oxt_xposdata_control(idx_active_control_oxt), oxt_yposdata_control(idx_active_control_oxt), 50, oxttreatment_color, 'filled');
        set(gca, 'Ydir', 'reverse');
        axis([-10 75 0 160])
        box off;
        axis off;
        
        print(sprintf('%s_%s%s', exp_type, exp_date, '-spatial_oxt_fromallcells'), '-depsc2');

        %plot spatial others
        figure(5);
        clf;
        subplot(1,2,2);
        scatter(others_xposdata_treatment, others_yposdata_treatment,  50, [0.5 0.5 0.5], 'filled');
        hold on;
        scatter(others_xposdata_treatment(idx_active_treatment_others), others_yposdata_treatment(idx_active_treatment_others), 50, otherstreatment_color, 'filled');
        set(gca, 'Ydir', 'reverse');
        axis([-10 75 0 160])
        box off;
        axis off;
        
        subplot(1,2,1);
        scatter(others_xposdata_control, others_yposdata_control,  50, [0.5 0.5 0.5], 'filled');
        hold on;
        scatter(others_xposdata_control(idx_active_control_others), others_yposdata_control(idx_active_control_others), 50, otherstreatment_color, 'filled');
        set(gca, 'Ydir', 'reverse');
        axis([-10 75 0 160])
        box off;
        axis off;
        
        print(sprintf('%s_%s%s', exp_type, exp_date,'-spatial_others_fromallcells'), '-depsc2');
        
    end
    
    %clear everything unnecessary
    clearvars -except alldata analyzeddata folderpath folders fnames_exp
    
end

save('alldissectedpERKdata.mat', 'analyzeddata', 'alldata')

%% NOW POOL DATA ACROSS ALL OPTIONS
close all;
f = 18;
oxtcontrol_color = 'k'
oxttreatment_color = [0 204/255 0] %[204/255 0 0];
otherscontrol_color = [0.5 0.5 0.5];
otherstreatment_color =[ 153/255 0 76/255] %[0 0 102/255];

nplot = 200; %number of neurons to downsample to. Can be modified

% if running from this section, may need to restate savefolder below
%change folder path accordingly
savefolder = '/Users/carolinewee/Desktop/pERK/epsfiles/'
cd(savefolder);

dataarray = struct;

fnames_exp = fieldnames(analyzeddata);

for e = [1:8,length(fnames_exp)]
    
    exp_type = char(fnames_exp(e));
    
    dataarray = setfield(dataarray, exp_type, []);
    
    fnames = fieldnames(alldata.(exp_type).data_control);
    
    pooled_oxt_controldata = [];
    pooled_oxt_treatmentdata = [];
    
    pooled_others_controldata = [];
    pooled_others_treatmentdata = [];
    
    pooled_oxt_controldata_FLUOR = [];
    pooled_oxt_treatmentdata_FLUOR = [];
    
    %ratios
    %oxt neurons using GFP
    pooled_ratios_control_FLUOR = [];
    pooled_ratios_treatment_FLUOR= [];
    
    %oxt neurons using all cells
    pooled_ratios_control_oxt = [];
    pooled_ratios_treatment_oxt = [];
    
    %other neurons using all cells
    pooled_ratios_control_others = [];
    pooled_ratios_treatment_others = [];
    
    %xyz
    pooled_oxt_xposdata_treatment = [];
    pooled_oxt_yposdata_treatment = [];
    pooled_oxt_xposdata_control = [];
    pooled_oxt_yposdata_control = []
    
    pooled_others_xposdata_control = [];
    pooled_others_yposdata_control = [];
    pooled_others_xposdata_treatment = [];
    pooled_others_yposdata_treatment = [];
    
    pooled_xposdata_treatment_FLUOR = [];
    pooled_yposdata_treatment_FLUOR = [];
    pooled_xposdata_control_FLUOR = [];
    pooled_yposdata_control_FLUOR = []
    
    for d = 1:length(fnames)
        
        exp_date = char(fnames(d));
        
        %pERK/tERK ratios (already normalized)
        pooled_oxt_controldata(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_controldata)) = analyzeddata.(exp_type).(exp_date).oxt_controldata;
        pooled_oxt_treatmentdata(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_treatmentdata)) = analyzeddata.(exp_type).(exp_date).oxt_treatmentdata;
        pooled_others_controldata(end+1:end+length(analyzeddata.(exp_type).(exp_date).others_controldata)) = analyzeddata.(exp_type).(exp_date).others_controldata;
        pooled_others_treatmentdata(end+1:end+length(analyzeddata.(exp_type).(exp_date).others_treatmentdata)) = analyzeddata.(exp_type).(exp_date).others_treatmentdata;
        pooled_oxt_controldata_FLUOR(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_controldata_FLUOR)) = analyzeddata.(exp_type).(exp_date).oxt_controldata_FLUOR;
        pooled_oxt_treatmentdata_FLUOR(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_treatmentdata_FLUOR)) = analyzeddata.(exp_type).(exp_date).oxt_treatmentdata_FLUOR;
        
        %x and y positions
        
        %oxt from GFP
        pooled_xposdata_control_FLUOR(end+1:end+length(analyzeddata.(exp_type).(exp_date).xposdata_control_FLUOR)) = analyzeddata.(exp_type).(exp_date).xposdata_control_FLUOR;
        pooled_yposdata_control_FLUOR(end+1:end+length(analyzeddata.(exp_type).(exp_date).yposdata_control_FLUOR)) = analyzeddata.(exp_type).(exp_date).yposdata_control_FLUOR;
        pooled_xposdata_treatment_FLUOR(end+1:end+length(analyzeddata.(exp_type).(exp_date).xposdata_treatment_FLUOR)) = analyzeddata.(exp_type).(exp_date).xposdata_treatment_FLUOR;
        pooled_yposdata_treatment_FLUOR(end+1:end+length(analyzeddata.(exp_type).(exp_date).yposdata_treatment_FLUOR)) = analyzeddata.(exp_type).(exp_date).yposdata_treatment_FLUOR;
        
        %oxt from all cells
        pooled_oxt_xposdata_control(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_xposdata_control)) = analyzeddata.(exp_type).(exp_date).oxt_xposdata_control;
        pooled_oxt_yposdata_control(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_yposdata_control)) = analyzeddata.(exp_type).(exp_date).oxt_yposdata_control;
        pooled_oxt_xposdata_treatment(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_xposdata_treatment)) = analyzeddata.(exp_type).(exp_date).oxt_xposdata_treatment;
        pooled_oxt_yposdata_treatment(end+1:end+length(analyzeddata.(exp_type).(exp_date).oxt_yposdata_treatment)) = analyzeddata.(exp_type).(exp_date).oxt_yposdata_treatment;
        
        %others from all cells
        pooled_others_xposdata_control(end+1:end+length(analyzeddata.(exp_type).(exp_date).others_xposdata_control)) = analyzeddata.(exp_type).(exp_date).others_xposdata_control;
        pooled_others_yposdata_control(end+1:end+length(analyzeddata.(exp_type).(exp_date).others_yposdata_control)) = analyzeddata.(exp_type).(exp_date).others_yposdata_control;
        pooled_others_xposdata_treatment(end+1:end+length(analyzeddata.(exp_type).(exp_date).others_xposdata_treatment)) = analyzeddata.(exp_type).(exp_date).others_xposdata_treatment;
        pooled_others_yposdata_treatment(end+1:end+length(analyzeddata.(exp_type).(exp_date).others_yposdata_treatment)) = analyzeddata.(exp_type).(exp_date).others_yposdata_treatment;
        
        %NORMALIZED pERK/tERK ratios by fish
        norm_ratios_treatment_FLUOR = analyzeddata.(exp_type).(exp_date).ratios_treatment_FLUOR;
        norm_ratios_control_FLUOR = analyzeddata.(exp_type).(exp_date).ratios_control_FLUOR;
        pooled_ratios_treatment_FLUOR = horzcat(pooled_ratios_treatment_FLUOR, norm_ratios_treatment_FLUOR);
        pooled_ratios_control_FLUOR = horzcat(pooled_ratios_control_FLUOR, norm_ratios_control_FLUOR);
        
        
        norm_ratios_treatment_oxt = analyzeddata.(exp_type).(exp_date).ratios_treatment_oxt; %divide by mean of all extracted neurons
        norm_ratios_control_oxt = analyzeddata.(exp_type).(exp_date).ratios_control_oxt;
        norm_ratios_treatment_others = analyzeddata.(exp_type).(exp_date).ratios_treatment_others; %divide by mean of all extracted neurons
        norm_ratios_control_others = analyzeddata.(exp_type).(exp_date).ratios_control_others;
        
        pooled_ratios_treatment_oxt = horzcat(pooled_ratios_treatment_oxt, norm_ratios_treatment_oxt);
        pooled_ratios_control_oxt = horzcat(pooled_ratios_control_oxt, norm_ratios_control_oxt);
        pooled_ratios_treatment_others = horzcat(pooled_ratios_treatment_others, norm_ratios_treatment_others);
        pooled_ratios_control_others = horzcat(pooled_ratios_control_others, norm_ratios_control_others);
        
    end
    
    % plot figure
    bandwidth = 0.15;
    ymax = 2
    figure(1);
    clf;
    
    [f_cont,xi_cont] = ksdensity(pooled_oxt_controldata, 'Bandwidth', bandwidth);
    [f_treatment,xi] = ksdensity(pooled_oxt_treatmentdata, 'Bandwidth', bandwidth);
    
    plot(xi_cont, f_cont, 'Color', oxtcontrol_color, 'LineWidth', 1);
    hold on;
    plot(xi, f_treatment, 'Color', oxttreatment_color, 'LineWidth', 1);
    set(gca, 'FontSize', f)
    box on;
    set(gca, 'FontSize', f)
    set(gca,'YTick',[])
    set(gca,'YTicklabel','')
    box off;
    hold on;
    
    [f_cont,xi_cont] = ksdensity(pooled_others_controldata, 'Bandwidth', bandwidth);
    [f_treatment,xi] = ksdensity(pooled_others_treatmentdata, 'Bandwidth', bandwidth);
    
    plot(xi_cont, f_cont, 'Color', otherscontrol_color, 'LineWidth', 1);
    hold on;
    plot(xi, f_treatment, 'Color', otherstreatment_color, 'LineWidth', 1);
    set(gca, 'FontSize', f)
    box on;
    set(gca, 'FontSize', f)
    set(gca,'YTick',[])
    set(gca,'YTicklabel','')
    %set(gca,'XTick',[])
    %set(gca,'XTicklabel','')
    axis([0 2.5 0 ymax])
    box off;
    
    print(sprintf('%s%s', exp_type, '-allcellshistogram'), '-depsc2');
    
    
    figure(2);
    clf;
    [f_cont,xi_cont] = ksdensity(pooled_oxt_controldata_FLUOR, 'Bandwidth', bandwidth);
    [f_treatment,xi] = ksdensity(pooled_oxt_treatmentdata_FLUOR, 'Bandwidth', bandwidth);
    
    % if using pdf rather than kernel density
    % f_cont = pdf(contdata, xi_cont);
    % f_treatment = pdf(treatmentdata, xi);
    
    plot(xi_cont, f_cont, 'Color', oxtcontrol_color, 'LineWidth', 1);
    hold on;
    plot(xi, f_treatment, 'Color', oxttreatment_color, 'LineWidth', 1);
    set(gca, 'FontSize', f)
    box on;
    set(gca, 'FontSize', f)
    set(gca,'YTick',[])
    set(gca,'YTicklabel','')
    %set(gca,'XTick',[])
    %set(gca,'XTicklabel','')
    box off;
    hold on;
    
    % this is same data from previous section
    [f_cont,xi_cont] = ksdensity(pooled_others_controldata, 'Bandwidth', bandwidth);
    [f_treatment,xi] = ksdensity(pooled_others_treatmentdata, 'Bandwidth', bandwidth);
    
    % f_cont = pdf(contdata, xi_cont);
    % f_treatment = pdf(treatmentdata, xi);
    
    plot(xi_cont, f_cont, 'Color', otherscontrol_color, 'LineWidth', 1);
    hold on;
    plot(xi, f_treatment, 'Color', otherstreatment_color, 'LineWidth', 1);
    set(gca, 'FontSize', f)
    box on;
    set(gca, 'FontSize', f)
    set(gca,'YTick',[])
    set(gca,'YTicklabel','')
    %set(gca,'XTick',[])
    %set(gca,'XTicklabel','')
    axis([0 2.5 0 ymax])
    box off;
    
    print(sprintf('%s%s', exp_type, '-allcellsandGFPhistogram'), '-depsc2');
    
    [h_ks, p_ks] = kstest2(pooled_oxt_controldata, pooled_oxt_treatmentdata) % oxt based on all cell data
    [p_rs, h_rs] = ranksum(pooled_oxt_controldata, pooled_oxt_treatmentdata)
    [h_ks2, p_ks2] = kstest2(pooled_others_controldata, pooled_others_treatmentdata)
    [p_rs2, h_rs2] = ranksum(pooled_others_controldata, pooled_others_treatmentdata)
    [h_ks3, p_ks3] = kstest2(pooled_oxt_controldata_FLUOR, pooled_oxt_treatmentdata_FLUOR)
    [p_rs3, h_rs3] = ranksum(pooled_oxt_controldata_FLUOR, pooled_oxt_treatmentdata_FLUOR)
    
    % also compare oxt treatment to others treatment
    [h_ks7, p_ks7] = kstest2(pooled_oxt_treatmentdata, pooled_others_treatmentdata);
    [p_rs7, h_rs7] = ranksum(pooled_oxt_treatmentdata, pooled_others_treatmentdata);
    [h_ks8, p_ks8] = kstest2(pooled_oxt_treatmentdata_FLUOR, pooled_others_treatmentdata);
    [p_rs8, h_rs8] = ranksum(pooled_oxt_treatmentdata_FLUOR, pooled_others_treatmentdata);
    
    % store statistics
    dataarray.(exp_type).stats.neurons.oxtfromallcells.kstest = [h_ks p_ks];
    dataarray.(exp_type).stats.neurons.oxtfromallcells.ranksum = [h_rs  p_rs];
    
    dataarray.(exp_type).stats.neurons.othersfromallcells.kstest = [h_ks2 p_ks2];
    dataarray.(exp_type).stats.neurons.othersfromallcells.ranksum = [h_rs2  p_rs2];
    
    dataarray.(exp_type).stats.neurons.oxtfromGFP.kstest = [h_ks3 p_ks3];
    dataarray.(exp_type).stats.neurons.oxtfromGFP.ranksum = [h_rs3  p_rs3];
    
    dataarray.(exp_type).stats.neurons.oxtvsothers.kstest = [h_ks7 p_ks7];
    dataarray.(exp_type).stats.neurons.oxtvsothers.ranksum = [h_rs7 p_rs7];
    dataarray.(exp_type).stats.neurons.oxtfromGFPvsothers.kstest = [h_ks8 p_ks8];
    dataarray.(exp_type).stats.neurons.oxtfromGFPvsothers.ranksum = [h_rs8 p_rs8];
    
    % store in data array
    dataarray.(exp_type).pooled_oxt_controldata = pooled_oxt_controldata;
    dataarray.(exp_type).pooled_oxt_treatmentdata = pooled_oxt_treatmentdata;
    dataarray.(exp_type).pooled_others_controldata = pooled_others_controldata;
    dataarray.(exp_type).pooled_others_treatmentdata = pooled_others_treatmentdata;
    dataarray.(exp_type).pooled_oxt_controldata_FLUOR = pooled_oxt_controldata_FLUOR;
    dataarray.(exp_type).pooled_oxt_treatmentdata_FLUOR = pooled_oxt_treatmentdata_FLUOR;
    dataarray.(exp_type).pooled_xposdata_treatment_FLUOR = pooled_xposdata_treatment_FLUOR;
    dataarray.(exp_type).pooled_yposdata_treatment_FLUOR = pooled_yposdata_treatment_FLUOR;
    dataarray.(exp_type).pooled_xposdata_control_FLUOR = pooled_xposdata_control_FLUOR;
    dataarray.(exp_type).pooled_yposdata_control_FLUOR = pooled_yposdata_control_FLUOR;
    
    dataarray.(exp_type).pooled_ratios_treatment_FLUOR = pooled_ratios_treatment_FLUOR;
    dataarray.(exp_type).pooled_ratios_control_FLUOR = pooled_ratios_control_FLUOR;
    dataarray.(exp_type).pooled_ratios_treatment_oxt = pooled_ratios_treatment_oxt;
    dataarray.(exp_type).pooled_ratios_control_oxt = pooled_ratios_control_oxt;
    dataarray.(exp_type).pooled_ratios_treatment_others = pooled_ratios_treatment_others;
    dataarray.(exp_type).pooled_ratios_control_others = pooled_ratios_control_others;
    
    % fraction responsive (normalized to mean of RESPECTIVE controls)
    pooled_all_controldata = horzcat(pooled_oxt_controldata,pooled_others_controldata);
    
    % find threshold in which only 3% of all control neurons would surpass X %
    % Uses findthreshold.m
    percentile = 10;
    active_threshold_pooled = findthreshold(pooled_all_controldata, percentile, 'above');
    suppressed_threshold_pooled = findthreshold(pooled_all_controldata, percentile, 'below');
    
    active_threshold_pooled_FLUOR = findthreshold(pooled_oxt_controldata_FLUOR, percentile, 'above');
    suppressed_threshold_pooled_FLUOR = findthreshold(pooled_oxt_controldata_FLUOR, percentile, 'below');
    
    % activated
    idx_active_control_others_pooled = find(pooled_others_controldata > active_threshold_pooled_FLUOR);
    idx_active_treatment_others_pooled = find(pooled_others_treatmentdata > active_threshold_pooled_FLUOR);
    idx_active_control_oxt_pooled = find(pooled_oxt_controldata > active_threshold_pooled_FLUOR);
    idx_active_treatment_oxt_pooled = find(pooled_oxt_treatmentdata > active_threshold_pooled_FLUOR);
    
    % inhibited
    idx_suppressed_control_others_pooled = find(pooled_others_controldata < suppressed_threshold_pooled_FLUOR);
    idx_suppressed_treatment_others_pooled = find(pooled_others_treatmentdata < suppressed_threshold_pooled_FLUOR);
    idx_suppressed_control_oxt_pooled = find(pooled_oxt_controldata < suppressed_threshold_pooled_FLUOR);
    idx_suppressed_treatment_oxt_pooled = find(pooled_oxt_treatmentdata < suppressed_threshold_pooled_FLUOR);
    
    fraction_responsive_others_control_pooled = length(idx_active_control_others_pooled)/length(pooled_others_controldata);
    fraction_responsive_oxt_control_pooled = length(idx_active_control_oxt_pooled)/length(pooled_oxt_controldata);
    fraction_responsive_others_treatment_pooled = length(idx_active_treatment_others_pooled)/length(pooled_others_treatmentdata);
    fraction_responsive_oxt_treatment_pooled = length(idx_active_treatment_oxt_pooled)/length(pooled_oxt_treatmentdata);
    
    fraction_suppressed_others_control_pooled = length(idx_suppressed_control_others_pooled)/length(pooled_others_controldata);
    fraction_suppressed_oxt_control_pooled = length(idx_suppressed_control_oxt_pooled)/length(pooled_oxt_controldata);
    fraction_suppressed_others_treatment_pooled = length(idx_suppressed_treatment_others_pooled)/length(pooled_others_treatmentdata);
    fraction_suppressed_oxt_treatment_pooled = length(idx_suppressed_treatment_oxt_pooled)/length(pooled_oxt_treatmentdata);
    
    dataarray.(exp_type).active_threshold_pooled= active_threshold_pooled;
    dataarray.(exp_type).fraction_responsive_others_control_pooled = fraction_responsive_others_control_pooled;
    dataarray.(exp_type).fraction_responsive_oxt_control_pooled = fraction_responsive_oxt_control_pooled;
    dataarray.(exp_type).fraction_responsive_others_treatment_pooled = fraction_responsive_others_treatment_pooled;
    dataarray.(exp_type).fraction_responsive_oxt_treatment_pooled = fraction_responsive_oxt_treatment_pooled;
    
    dataarray.(exp_type).suppressed_threshold_pooled= suppressed_threshold_pooled;
    dataarray.(exp_type).fraction_suppressed_others_control_pooled = fraction_suppressed_others_control_pooled;
    dataarray.(exp_type).fraction_suppressed_oxt_control_pooled = fraction_suppressed_oxt_control_pooled;
    dataarray.(exp_type).fraction_suppressed_others_treatment_pooled = fraction_suppressed_others_treatment_pooled;
    dataarray.(exp_type).fraction_suppressed_oxt_treatment_pooled = fraction_suppressed_oxt_treatment_pooled;
    
    % fraction responsive OXT using GFP fluorescence
    active_threshold_pooled_FLUOR = findthreshold(pooled_oxt_controldata_FLUOR, percentile, 'above');
    suppressed_threshold_pooled_FLUOR = findthreshold(pooled_oxt_controldata_FLUOR, percentile, 'below');
    
    idx_active_control_FLUOR_pooled = find(pooled_oxt_controldata_FLUOR>active_threshold_pooled_FLUOR);
    idx_active_treatment_FLUOR_pooled = find(pooled_oxt_treatmentdata_FLUOR>active_threshold_pooled_FLUOR);
    
    idx_suppressed_control_FLUOR_pooled = find(pooled_oxt_controldata_FLUOR < suppressed_threshold_pooled_FLUOR);
    idx_suppressed_treatment_FLUOR_pooled = find(pooled_oxt_treatmentdata_FLUOR < suppressed_threshold_pooled_FLUOR);
    
    fraction_responsive_oxt_control_FLUOR_pooled = length(idx_active_control_FLUOR_pooled)/length(pooled_oxt_controldata_FLUOR);
    fraction_responsive_oxt_treatment_FLUOR_pooled = length(idx_active_treatment_FLUOR_pooled)/length(pooled_oxt_treatmentdata_FLUOR);
    dataarray.(exp_type).active_threshold_pooled_FLUOR = active_threshold_pooled_FLUOR;
    dataarray.(exp_type).fraction_responsive_oxt_control_FLUOR_pooled = fraction_responsive_oxt_control_FLUOR_pooled;
    dataarray.(exp_type).fraction_responsive_oxt_treatment_FLUOR_pooled = fraction_responsive_oxt_treatment_FLUOR_pooled;
    
    fraction_suppressed_oxt_control_FLUOR_pooled = length(idx_suppressed_control_FLUOR_pooled)/length(pooled_oxt_controldata_FLUOR);
    fraction_suppressed_oxt_treatment_FLUOR_pooled = length(idx_suppressed_treatment_FLUOR_pooled)/length(pooled_oxt_treatmentdata_FLUOR);
    dataarray.(exp_type).suppressed_threshold_pooled_FLUOR = suppressed_threshold_pooled_FLUOR;
    dataarray.(exp_type).fraction_suppressed_oxt_control_FLUOR_pooled = fraction_suppressed_oxt_control_FLUOR_pooled;
    dataarray.(exp_type).fraction_suppressed_oxt_treatment_FLUOR_pooled = fraction_suppressed_oxt_treatment_FLUOR_pooled;
    
    % PLOT SPATIAL.
    % FIRST RESAMPLE SO THAT sAME NUMBER OF NEURONS FOR ALL STIMULI
    % Use function resampleindices.m
    
    % plot spatial oxt fluor
    figure(3);
    clf;
    [idx_subset, new_idx] = resampleindices(pooled_oxt_treatmentdata_FLUOR, idx_active_treatment_FLUOR_pooled, nplot);
    
    subplot(2,4,2); %treatment
    scatter(pooled_xposdata_treatment_FLUOR(idx_subset), pooled_yposdata_treatment_FLUOR(idx_subset),  50, [0.5 0.5 0.5], 'filled');
    hold on;
    scatter(pooled_xposdata_treatment_FLUOR(new_idx), pooled_yposdata_treatment_FLUOR(new_idx), 50, oxttreatment_color, 'filled');
    %scatter(pooled_xposdata_treatment_FLUOR(idx_suppressed_treatment_FLUOR_pooled), pooled_yposdata_treatment_FLUOR(idx_suppressed_treatment_FLUOR_pooled), 50, 'k', 'filled');
    set(gca, 'Ydir', 'reverse');
    axis([-10 75 0 160])
    box off;
    axis off;
    
    [idx_subset, new_idx] = resampleindices(pooled_oxt_controldata_FLUOR, idx_active_control_FLUOR_pooled, nplot);
    
    subplot(2,4,1); %control
    scatter(pooled_xposdata_control_FLUOR(idx_subset), pooled_yposdata_control_FLUOR(idx_subset),  50, [0.5 0.5 0.5], 'filled');
    hold on;
    scatter(pooled_xposdata_control_FLUOR(new_idx), pooled_yposdata_control_FLUOR(new_idx), 50, oxttreatment_color, 'filled');
    %scatter(pooled_xposdata_control_FLUOR(idx_suppressed_control_FLUOR_pooled), pooled_yposdata_control_FLUOR(idx_suppressed_control_FLUOR_pooled), 50, 'k', 'filled');
    set(gca, 'Ydir', 'reverse');
    axis([-10 75 0 160])
    box off;
    axis off;
    
    print(sprintf('%s-%d%s', exp_type, nplot, '-spatial_oxt_fromGFP'), '-depsc2');
    
    % plot spatial oxt from all cells
    
    figure(4);
    clf;
    
    [idx_subset, new_idx] = resampleindices(pooled_oxt_treatmentdata, idx_active_treatment_oxt_pooled, nplot);
    
    subplot(2,4,2); %treatment
    scatter(pooled_oxt_xposdata_treatment(idx_subset), pooled_oxt_yposdata_treatment(idx_subset),  50, [0.5 0.5 0.5], 'filled');
    hold on;
    scatter(pooled_oxt_xposdata_treatment(new_idx), pooled_oxt_yposdata_treatment(new_idx), 50, oxttreatment_color, 'filled');
    %scatter(pooled_oxt_xposdata_treatment(idx_suppressed_treatment_oxt_pooled), pooled_oxt_yposdata_treatment(idx_suppressed_treatment_oxt_pooled), 50, 'k', 'filled');
    set(gca, 'Ydir', 'reverse');
    axis([-10 75 0 160])
    box off;
    axis off
    
    [idx_subset, new_idx] = resampleindices(pooled_oxt_controldata, idx_active_control_oxt_pooled, nplot);
    
    subplot(2,4,1);
    scatter(pooled_oxt_xposdata_control(idx_subset), pooled_oxt_yposdata_control(idx_subset),  50, [0.5 0.5 0.5], 'filled');
    hold on;
    scatter(pooled_oxt_xposdata_control(new_idx), pooled_oxt_yposdata_control(new_idx), 50, oxttreatment_color, 'filled');
    %scatter(pooled_oxt_xposdata_control(idx_suppressed_control_oxt_pooled), pooled_oxt_yposdata_control(idx_suppressed_control_oxt_pooled), 50, 'k', 'filled');
    set(gca, 'Ydir', 'reverse');
    axis([-10 75 0 160])
    box off;
    axis off;
    
    print(sprintf('%s-%d%s', exp_type, nplot, '-spatial_oxt_fromallcells'), '-depsc2');
    
    % other neurons
    
    nplot_others = 800;
    
    figure(5);
    clf;
    
    [idx_subset, new_idx] = resampleindices(pooled_others_treatmentdata, idx_active_treatment_others_pooled,  nplot_others);
    
    subplot(2,4,2);
    scatter(pooled_others_xposdata_treatment(idx_subset), pooled_others_yposdata_treatment(idx_subset),  50, [0.5 0.5 0.5], 'filled');
    hold on;
    scatter(pooled_others_xposdata_treatment(new_idx), pooled_others_yposdata_treatment(new_idx), 50, otherstreatment_color, 'filled');
    % scatter(pooled_others_xposdata_treatment(idx_suppressed_treatment_others_pooled), pooled_others_yposdata_treatment(idx_suppressed_treatment_others_pooled), 50, 'k', 'filled');
    set(gca, 'Ydir', 'reverse');
    axis([-10 75 0 160])
    box off;
    axis off;
    
    [idx_subset, new_idx] = resampleindices(pooled_others_controldata, idx_active_control_others_pooled, nplot_others);
    
    subplot(2,4,1);
    scatter(pooled_others_xposdata_control(idx_subset), pooled_others_yposdata_control(idx_subset),  50, [0.5 0.5 0.5], 'filled');
    hold on;
    scatter(pooled_others_xposdata_control(new_idx), pooled_others_yposdata_control(new_idx), 50, otherstreatment_color, 'filled');
    % scatter(pooled_others_xposdata_control(idx_suppressed_control_others_pooled), pooled_others_yposdata_control(idx_suppressed_control_others_pooled), 50, 'k', 'filled');
    set(gca, 'Ydir', 'reverse');
    axis([-10 75 0 160])
    box off;
    axis off;
    
    print(sprintf('%s-%d%s', exp_type, nplot_others, '-spatial_others_fromallcells'), '-depsc2');
    
    %% scatter plot by fish
    
    % from all cells
    pooled_ratios_control_oxt_byfish = nanmean(pooled_ratios_control_oxt,1);
    pooled_ratios_treatment_oxt_byfish = nanmean(pooled_ratios_treatment_oxt,1);
    pooled_ratios_control_others_byfish = nanmean(pooled_ratios_control_others,1);
    pooled_ratios_treatment_others_byfish = nanmean(pooled_ratios_treatment_others,1);
    
    % regardless of neuron types
    nfish_control = length(pooled_ratios_control_oxt_byfish(~isnan(pooled_ratios_control_oxt_byfish)));
    nfish_treatment = length(pooled_ratios_treatment_oxt_byfish(~isnan(pooled_ratios_treatment_oxt_byfish)));
    
    [p h] = ranksum(pooled_ratios_control_oxt_byfish, pooled_ratios_treatment_oxt_byfish);
    [p2 h2] = ranksum(pooled_ratios_control_others_byfish, pooled_ratios_treatment_others_byfish);
    
    dataarray.(exp_type). pooled_ratios_control_oxt_byfish = pooled_ratios_control_oxt_byfish;
    dataarray.(exp_type). pooled_ratios_treatment_oxt_byfish = pooled_ratios_treatment_oxt_byfish;
    dataarray.(exp_type). pooled_ratios_control_others_byfish = pooled_ratios_control_others_byfish;
    dataarray.(exp_type). pooled_ratios_treatment_others_byfish = pooled_ratios_treatment_others_byfish;
    dataarray.(exp_type).stats.fish.oxtfromallcells.ranksum = [p h nanmean(pooled_ratios_control_oxt_byfish) nanmean(pooled_ratios_treatment_oxt_byfish) nfish_control nfish_treatment];
    dataarray.(exp_type).stats.fish.others.ranksum = [p2 h2 nanmean(pooled_ratios_control_others_byfish) nanmean(pooled_ratios_treatment_others_byfish) nfish_control nfish_treatment];
    
    %from GFP data
    pooled_ratios_control_byfish = nanmean(pooled_ratios_control_FLUOR,1);
    pooled_ratios_treatment_byfish = nanmean(pooled_ratios_treatment_FLUOR, 1);
    
    [p3 h3] = ranksum(pooled_ratios_control_byfish, pooled_ratios_treatment_byfish);
    
    dataarray.(exp_type). pooled_ratios_control_byfish = pooled_ratios_control_byfish;
    dataarray.(exp_type). pooled_ratios_treatment_byfish = pooled_ratios_treatment_byfish;
    dataarray.(exp_type).stats.fish.oxtfromGFP.ranksum = [p3 h3 nanmean(pooled_ratios_control_byfish) nanmean(pooled_ratios_treatment_byfish) nfish_control nfish_treatment];
    
    lightgray = [0.5 0.5 0.5];
    darkgray = [0.2 0.2 0.2];
    marker_size = 100;
    figure(6)
    clf;
    
    % normalize to mean of control
    ratios_others =  pooled_ratios_treatment_others_byfish/nanmean(pooled_ratios_control_others_byfish);
    ratios_oxt = pooled_ratios_treatment_byfish/nanmean(pooled_ratios_control_byfish);
    
    % to plot means
    bar([1 2], [nanmean(ratios_oxt) nanmean(ratios_others)], 'k');
    hold on;
    scatter(ones(1, length(pooled_ratios_treatment_byfish)), ratios_oxt , marker_size, oxttreatment_color, 'filled');
    scatter(2*ones(1, length(pooled_ratios_treatment_others_byfish)), ratios_others , marker_size, otherstreatment_color, 'filled');
    box off;
    axis([0.5 2.5 0.8 1.7]);
    
    %     scatter(ones(1, length(pooled_ratios_control_others_byfish)), pooled_ratios_control_others_byfish, marker_size, otherscontrol_color, 'filled');
    %     hold on;
    %     scatter(2*ones(1, length(pooled_ratios_treatment_others_byfish)), pooled_ratios_treatment_others_byfish, marker_size, otherstreatment_color, 'filled');
    %
    %     scatter(1.05*ones(1, length(pooled_ratios_control_byfish)), pooled_ratios_control_byfish, marker_size, oxtcontrol_color, 'filled');
    %     scatter(2.05*ones(1, length(pooled_ratios_treatment_byfish)), pooled_ratios_treatment_byfish, marker_size, oxttreatment_color, 'filled');
    
    errorbar([1], nanmean(ratios_oxt), nanstd(ratios_oxt)/sqrt(nfish_treatment), 'color', 'k');
    errorbar([2], nanmean(ratios_others), nanstd(ratios_others)/sqrt(nfish_treatment), 'color',  'k');
    
    print(sprintf('%s%s', exp_type, '-scatterbyfish'), '-depsc2');
    
    
end

cd(savefolder);
save('alldissectedpERKdata.mat', 'analyzeddata', 'alldata', 'dataarray')

