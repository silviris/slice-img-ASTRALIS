# slice-img-ASTRALIS

Analysis pipeline
0) prepare folder/save info:
access 'C:\Users\silvi\Desktop\troubleShoot' and copy its content in the slice folder
place raw images folders in 'raw'
open '.xml' file to have pxSize and frRate at hand
save all relevant mouse/slice/expt info using matlab GUI
1) ImageJ:
extract videos from folders, calibrate (pxSize), adjust brightness, and save to 'videos' folder [use macro:'loadCalibSave_SP.ijm']
extract full fov fluorescence over time, find maxima, draw circular ROIs, and extract traces to separate folders in '\ANALYSIS_ImageJ\'>  -> '\ANALYSIS_imageJ\fovFluo' & '..\somataFluo'
use Matlab scripts to extract fov fluorescence and plot it for all videos, as well as saving to the slice folder the '\results.mat' file, and to the 'ImageJ\fovFluo' folder the traces and peak detection figures
2) AQuA:
copy videos to analyze to '\ANALYSIS_aqua\'
run aqua on one of the videos (ideally the stim video) to find the optimal parameters for that slice, and save them as 'parameters1.csv' in '\ANALYSIS_aqua\'
save masks or landmarks as 'cellMask.mat/lmkMask.mat' in '\ANALYSIS_aqua\'
perform batch detection, and then batch fix_events

NB: from ImageJ and Aqua analysis -> output "results.mat" file in sliceDir, containing:
- R structure: info & results
- R_cnd: order of conditions by video number

3) T table ('allSlicesSummary.m'):
collect all data in R and R_cnd to create a table for plotting data from all slices easily
4) plotting from T:

5) pooling data across slices:


0) prepare folder/save info:
% rootDir = 'C:\Users\silvi\Documents\MATLAB';
% % you can open the script containing these lines here>
% open('C:\Users\silvi\Documents\MATLAB\silvia\templates_SP\userInputINFO_mouseSliceExpt.mlx')
clear
close all
 input MOUSE INFO (save in mouse folder)
run('C:\Users\silvi\Documents\MATLAB\silvia\GUIinfoMouse\mouseInfoGUI.m')
 input IMAGING PARAMETERS (saved in SLICE folder)  [NB: to save time, duplicate to other slice folders if img parameters are unchanged! ]
run('C:\Users\silvi\Documents\MATLAB\silvia\GUIinfoMouse\imgInfoGUI.m')
input SLICE INFO (saved in slice folder)
run('C:\Users\silvi\Documents\MATLAB\silvia\GUIinfoMouse\sliceInfoGUI.m')
input VIDEO INFO (saved in slice folder, with index indicating video number [eg: '001_videoInfo.mat'])
NB: you can use save and add button for the next video of this slice
run('C:\Users\silvi\Documents\MATLAB\silvia\GUIinfoMouse\videoInfoGUI.m')
NNB: for the next slice, you will want to remove 'sliceInfo' from the base workspace; same for the next mouse, remove 'mouseInfo' from workspace
% clear sliceInfo 
% clear mouseInfo
Create 'R' & 'R_cnd'
update cndID in 'xxx_videoInfo'
clearvars -except mouseInfo imgInfo sliceInfo
% update videoInfo with looped cndID>
create_CNDloop(char(sliceInfo.sliceDir))
create structures R and R_cnd and save them in "results.mat" in sliceDir>
NB: you can run this with no var in the ws and then a popup will ask you to select the mouse and then the slice folder
% collect info and save in "results.mat" in sliceDir
storeR_info(mouseInfo,sliceInfo,imgInfo)

1) ImageJ:
After using the ImageJ macro 'loadCalibSave_SP.ijm', export the traces and plot them with the script 'templates_SP\imageJdata_findPeakMaxOnsetOffset.mlx'
close all
clearvars -except mouseInfo imgInfo sliceInfo
Retrieve data:
sliceDir = char(sliceInfo.sliceDir);
dataDir = [sliceDir,'\ANALYSIS_ImageJ\fovFluo'];
% retrieve pxSize and frRate from 'imgInfo.m'>
frRate = imgInfo.frRate;
pxSize = imgInfo.pxSize;
Set parameters for detection of peaks in 'full-fov' data:
% set window (in seconds) for calculating dF/F>
dFF_window = 30; % seconds
% set by which CND you want to calculate the "noise" of all traces>
norm_CND = 'bl';
% set threshold for peak detection> {will take whichever is greater}
% how many std above baseline?
stdN= 3; 
% peaks must be at least over 5% change in dFF!
peakThreshold = 0.05; % 0.05
% correction factor for threshold> must be at least above x times the noise
% of the full baseline video>
timesNoise = 10;
minDist = 10; %s
% Smooth the data using the loess method with a span of 10% of
% the data>
spanSmooth = .1;
detect the peaks in all 'full-fov' traces for that slice:
run(['C:\Users\silvi\Documents\MATLAB\silvia\templates_SP\' ...
    'imageJdata_findPeakMaxOnsetOffset.mlx'])
save figures?
saveFigAsSvg(allTraces_f ,"traces_fov",dataDir,false)
saveFigAsSvg(peakDetect_f,"peakDetection_fov",dataDir,false)
savefig(peakDetect_f, [dataDir,'\peakDetection_fov']) % saved as matlab fig to allow navigation
saveFigAsSvg(peakQuantif_f,"respQuantif_fov",dataDir,false)

clear allTraces_f peakDetect_f peakQuantif_f
close all
store ImageJ data:
save([dataDir,'\resImageJ.mat'],'res') % in ImageJ dir
save([dataDir,'\detectionParam.mat'],'detectionParam') % in ImageJ dir
% % load( [sliceDir,'\ANALYSIS_ImageJ\fovFluo\resImageJ.mat']) % use to load res
storeR_fov(sliceDir,res) % ..as R, R_cnd in sliceDir ('results.mat')

clearvars -except mouseInfo imgInfo sliceInfo

2) AQuA:
2.1) aqua event detection:
aquaGuiFolder = 'C:\Users\silvi\Documents\MATLAB\aqua-19Dec19';
run aqua to find the optimal parameters for that slice:
aqua_gui
NB: I usually use a lower voxel Int Threshold than signal Int. Thresh. (by .1 or .2)
NB: if there are masks or landmarks, save them as> 'cellMask.mat' // 'lmkMask.mat' in '\ANALYSIS_aqua\'
after you found the parameters you want, save 'parameters1.csv' in '\ANALYSIS_aqua\' for batch detection. [Leave just that 1 column!]
process all the videos with the same parameters using the batch processing script, and save the output to '\ANALYSIS_aqua\':
run([aquaGuiFolder,'\aqua_cmd_batch.m'])
* OPTIONAL *: for experiments where the onset is really critical, fix events using "batch_fixEvents_SP.m"
run([aquaGuiFolder,'\batch_fixEvents_SP.m'])

2.2) spatio-temporal  plots, extract aqua features, save>
plot> 'silvia\plotsSP\spRasterFromPiaPropAP_DA.m'
save> 'R/R_cnd' in 'results.mat':
clearvars -except mouseInfo imgInfo sliceInfo
sliceDir = char(sliceInfo.sliceDir);
aquaFiles = dir(fullfile([sliceDir,'\ANALYSIS_aqua\'],'*_aqua.mat'));  % nPlots = length(aquaFiles);

for k=1:length(aquaFiles)
%     currentCnd = R_cnd{k};
    f0 = aquaFiles(k).name;
   
    % get cndID
    tmp = strsplit(f0,'_'); 
    cndIdx = tmp{3}; clear tmp 
 
    %load expt info>
    load([sliceDir,'\',cndIdx,'_videoInfo.mat']) 
    
    % load aqua file>
    load([sliceDir,'\ANALYSIS_aqua\',f0]); clear f0

extract aqua features>
    [aqua] = exportAquaRes(res, cndIdx);
    allCND.(cndID).aqua = aqua; % cndID is loaded from '_videoInfo.mat'
plot spatial raster>
    xRange = [0 600];
    run(['C:\Users\silvi\Documents\MATLAB\silvia\plotsSP\' ...
         'spRasterFromPiaPropAP_DA'])
    set(gcf,'Visible','on')
    sgtitle(cndID)
    saveFigAsSvg(f,[cndIdx,'_spTraster',],[sliceDir,'\ANALYSIS_aqua\'],false);
    
    clearvars -except aquaFiles allCND sliceDir k mouseInfo imgInfo sliceInfo
end
store extracted aqua data in [sliceDir,'\results.mat'] and in [aquaDir,'\allCND.mat']>
close all; clear aquaFiles k
storeR_aqua(sliceDir,allCND)   % save in 'results.m' in sliceDir
save([sliceDir,'\ANALYSIS_aqua\allCND.mat'],'allCND'); clear allCND % in Aqua dir

2.3)  filter events:
clearvars -except sliceInfo mouseInfo imgInfo sliceDir
load([sliceDir,'\results.mat']) % for R R_cnd
extract only aqua events that happen in a 3/5 min window after tON+perfDelay:
aquaWindowLength = [180 300]; % in seconds
find event whose onset time falls within the window from t=0 or t=drug, then save figures >
for ii=1:numel(aquaWindowLength)
    % NB: this script adds 'aqua.evFlt180s' & 'aqua.evFlt300s' to R, to
    % store the ID of events filtered using those windows>
    run(['C:\Users\silvi\Documents\MATLAB\silvia\tmp\' ...
         'plotAquaMean1slice.mlx'])
    % save figures and filtered analysis>
    saveFigAsSvg(aquaQuantif_f,['aquaVarMean_wind',num2str(aquaWindowLength(ii)),'s'],...
                 [sliceDir,'\ANALYSIS_aqua'],false)
    close(aquaQuantif_f)
end
store 'flt' filtered means of variables in [aquaDir,'\allCND_flt_means.mat'], where the rows of 'flt' correspond to the length of each  filtering windows, the fields of 'flt' correspond to event variables, and columns within fields are the different conditions (as in 'cnd') >
cnd = R_cnd;
save([sliceDir,'\ANALYSIS_aqua\'...
     'allCND_flt_evProp_means.mat'],... % filename
     'flt',...  main structure
     'aquaWindowLength',... length of each filtering windows (=rown in 'flt') [in sec]
     'plotVars',... (= fields of 'flt')
     'varLabels',... (= labels for 'plotVars')
     'cnd') % corresponding R_cnd for columns in each field of 'flt'
save back the IDs of the filtered events to 'fltEventsIDs' in the R variable by overriding '\results.m'>
save([sliceDir,'\results.mat'],'R','R_cnd')

clearvars -except sliceInfo mouseInfo imgInfo...
                  sliceDir ...
                  R R_cnd ...
                  aquaWindowLength

2.4) calculate general aqua results>
which CND is drug?
% drugIDX = find(contains(R_cnd, 'DA'));
what CND are you using for normalization?
% blIDX = find(strcmp(R_cnd, 'bl'));
which window length to use?
% w=1;  % [1==180s; 2==300s]
calculate the difference in response between drug and baseline, normalized by baseline, using these variables>
% deltaVars = {'evN','latency','dur','area'};
% run(['C:\Users\silvi\Documents\MATLAB\silvia\tmp\' ...
%     'deltaBlDrug_aqua'])
% 
save Δbl-drug in '\ANALYSIS_aqua'>
% save([sliceDir,'\ANALYSIS_aqua\delta_blDrug_AquaFlt',num2str(aquaWindowLength(w)),'s.mat'],...
%      '-struct','delta_blDrug')  
% % === the following done in next section> =====
% % % ... add aqua results to table 'T'>

2.5) calculate probActivePx for that window>
clearvars -except mouseInfo imgInfo sliceInfo sliceDir R R_cnd aquaWindowLength
normalize by how many times a px is active at any CND or within each video?
throughCND = true; % [true==through CND; false==within CND]
calculate 'pActivePx', normalize through CND, plot map/trace and trace-CDF, and save figures>
for w=1:numel(aquaWindowLength)
    run(['C:\Users\silvi\Documents\MATLAB\silvia\tmp\'...
    'probActivePx_acrossCND.mlx'])
    close all
save maps-traces 'pActivePx' in '\ANALYSIS_aqua'>
    save([sliceDir,'\ANALYSIS_aqua\pActivePx_tracesMaps_w',num2str(aquaWindowLength(w)),'s.mat'],...
         'pActive_trace','pActive_map') 
calculate mean pActivePx and save in 'allCND_flt_means.m' > {will be used later to populate T table}
    flt(w).pActPx = cellfun(@mean, pActive_trace);
    clear pActive_trace pActive_map
end

save([sliceDir,'\ANALYSIS_aqua\'...
     'allCND_flt_pActPx_means.mat'],...
     'flt') 

3) populate T table ('allSlicesSummary.m'):
clearvars -except sliceDir %sliceDir mouseInfo imgInfo aquaWindowLength w
% sliceDir = char(sliceInfo.sliceDir);

% ==========TABLE T===========
% get variable names from Table>
% var_names = T.Properties.VariableNames;
% save('C:\Users\silvi\Documents\MATLAB\silvia\templates_SP\varNames\T_varNames.mat','var_names')
% load('C:\Users\silvi\Documents\MATLAB\silvia\templates_SP\varNames\T_varNames.mat')
% % create empty T table with var names 'var_names'>
% var_empty = cell(1,length(var_names));
% T = cell2table(var_empty, 'VariableNames', var_names);
% save('allSlicesSummary','T')
% =============================
*** check also 'table_cellToMatNaNs.m' for tips on how to use tables ***
load R & R_cnd with all info & ImageJ results>
load([sliceDir,'\results.mat'])
load the summary of aqua data>
ev_means = load([sliceDir,'\ANALYSIS_aqua\allCND_flt_evProp_means.mat']);
pActPx_means = load([sliceDir,'\ANALYSIS_aqua\allCND_flt_pActPx_means.mat']);
load T>
tmp = load('allSlicesSummary.mat'); oldT = tmp.T; clear tmp 
% % use this if no table exists>  %oldT = [];
store the info/ImageJ/aqua data for this slice in the table 'T' from R & R_cnd:
sliceIdx = []; %[]
[T, sliceIdx] = addRow_allSlicesSummary(R, R_cnd, ev_means, pActPx_means, oldT,sliceIdx);
clear oldT ev_means pActPx_means
% % to add just mouse/slice/img info use this>
% [T, sliceIdx] = addRow_allSlicesSummary([],[],[],[],oldT,[]); clear oldT
save back the slice ID to the R variable>
for k = 1:length(R_cnd), R.(R_cnd{k}).sliceID = sliceIdx; end 
save([sliceDir,'\results.mat'],'R','R_cnd');
save table T>
T = sortrows(T); % sort by sliceID
save('allSlicesSummary','T')
clear T k sliceIdx R R_cnd

4) plotting from T:
% ..also for table T, then plot 'sex','age','coordinates'>>
open('plottingFromTables.mlx')

5) pooling data across slices:
% need to loop through certain slices/CND
open('C:\Users\silvi\Documents\MATLAB\silvia\tmp\storeAllSlicesFromT.m')
% check past code for doing this!! (Fluo4 analysis/IP3R2 ko)

% to plot all events across slices/cnd filtered in a certain window>
open('C:\Users\silvi\Documents\MATLAB\silvia\tmp\plotAllSlices_aquaEvents.m')

% to plot pActPx traces as imagesc>
open('C:\Users\silvi\Documents\MATLAB\silvia\tmp\plotAllSlices_pActTrace.m')

