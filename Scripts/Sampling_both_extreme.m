clear 
clc
Rf_Acc = 24;
compound_window = 1.5; % half / days
dec_tim = 1.5;          % half /days

path='H:\OneDrive - University of Central Florida\Analysis_AR_two_way_8_events_splitted_marg_comb_SS_Hmo\';
load('H:\OneDrive - University of Central Florida\Analysis_AR_two_way_8_events_splitted_marg_comb_SS_Hmo\Creating_WLONC\Time_series_data.mat');


%%
% RF_thres=nan(77,1);
% NTTWL_thres = nan(77,1);
% 
% for loc = 3%1:numel(shp)
%     Data = Time_series_data{loc,1};
% %     RF_acc_temp = movsum(Data(:,2), [Rf_Acc - 1, 0], 'omitnan');
% %     [~, v_peaks_RF] = decluster_POT(Data(:,1), RF_acc_temp, 1, 3, 1);
% %     RF_thres(loc,1) = prctile(v_peaks_RF,50);   
%     [~, v_peaks_NTTWL] = decluster_POT(Data(:,1), [Data(:,5)+Data(:,6)], 1, 3, 1);
%     NTTWL_thres(loc,1) = prctile(v_peaks_NTTWL,50);
% end

%%
RF_percentile = 90;
SS_percentile = 90;
POT_both_extreme={};
POT_RF_only_extreme={};
POT_NTTWL_only_extreme={};

SED_RF=[];
SED_WLAM=[];
SED_SS=[];
for loc = 1:numel(shp)
    Data = Time_series_data{loc,1};
    
    %RF_acc_temp = movsum(Data(:,2), [Rf_Acc - 1, 0], 'omitnan');
    %RF_acc_temp(RF_acc_temp==0,:)=[];
    %threshold_RF = prctile(RF_acc_temp,RF_percentile); 

    % new method: caculate the daily RF and get the 80th percentile
    dateDates = floor(Data(:,1));
    [uniqueDays, ~, idx] = unique(dateDates);
    dailyRainfall = [uniqueDays, accumarray(idx, Data(:,2))];

    threshold_RF = prctile(dailyRainfall(:,2),RF_percentile); 
    threshold_WL = prctile([Data(:,5)+Data(:,6)],SS_percentile);
    threshold_SS = prctile(Data(:,5),SS_percentile);
    
%     [SED_RF(loc,1),SED_RF(loc,2)]=find_parameters([Data(:,1),Data(:,2)],threshold_RF);
%     [SED_WLAM(loc,1),SED_WLAM(loc,2)]=find_parameters([Data(:,1),[Data(:,5)+Data(:,6)]],threshold_WL);
%     [SED_SS(loc,1),SED_SS(loc,2)]=find_parameters([Data(:,1),Data(:,5)],threshold_SS);

    
    out = find_POT_compound_RF_WLNC(Data, threshold_RF, threshold_WL, 'RfAccHours', Rf_Acc, 'compound_window', compound_window, 'dec_tim', dec_tim, 'doPlots', false);
    POT_both_extreme{loc,1}=out.POT;   % [t_RF, AccRF, t_Hm0, Hm0]
    POT_both_extreme{loc,2}=threshold_RF;
    POT_both_extreme{loc,3}=threshold_WL;

    out2 = find_POT_oneway_RF(Data, threshold_RF, 'RfAccHours', Rf_Acc, 'compound_window', compound_window, 'dec_tim', dec_tim, 'doPlots', false);

    out3 = find_POT_oneway_NTTWL(Data, threshold_WL, 'compound_window', compound_window, 'dec_tim', dec_tim, 'doPlots', false);

    POT_RF_only_extreme{loc,1}=out2.POT;   % [t_RF, AccRF]
    POT_RF_only_extreme{loc,2}=threshold_RF;
    POT_NTTWL_only_extreme{loc,1}=out3.POT;   % [t_RF, AccRF]
    POT_NTTWL_only_extreme{loc,2}=threshold_WL;



end
%%
save('POT_both_extreme.mat','POT_both_extreme','POT_RF_only_extreme',"POT_NTTWL_only_extreme");

close all
%%
load POT_both_extreme.mat
load CoastPieces.mat
load coastBuf.mat
% loading the coastline
coastline = shaperead("H:\OneDrive - University of Central Florida\Analysis_AR_two_way\Coastline\ne_10m_coastline\ne_10m_coastline.shp","BoundingBox",[-128, 22, ;-102, 49, ]);


  hFig = plotCoastPiecesColored(CoastPieces, SED_RF(:,1), [30 51], [-130 -110], ...
      'Title', 'SED of RF', ...
      'ColorbarLabel', 'SED (h)', ...
      'Colormap', turbo(256), ...
      'LineWidth', 4, ...
      'FontSize', 12, ...
      'BaseMap', 'grayland', ...
      'SaveFig', true, ...
      'SavePNG', true, ...
      'DPI', 300, ...
      'Coastline', coastline, ...
      'FigWidth', 8, ...
      'FigHeight', 7, ...
      'InsetBox1',[33.28 34.58 -118.5 -118], ...
      'InsetPos1', [0.6 0.25 0.2 0.3], ...
      'InsetBox2',[36.83 38.62 -123.16 -121.84], ...
      'InsetPos2', [0.6 0.6 0.2 0.3]);



  hFig = plotCoastPiecesColored(CoastPieces, SED_WLAM(:,1), [30 51], [-130 -110], ...
      'Title', 'SED of EWLAN', ...
      'ColorbarLabel', 'SED (h)', ...
      'Colormap', turbo(256), ...
      'LineWidth', 4, ...
      'FontSize', 12, ...
      'BaseMap', 'grayland', ...
      'SaveFig', true, ...
      'SavePNG', true, ...
      'DPI', 300, ...
      'Coastline', coastline, ...
      'FigWidth', 8, ...
      'FigHeight', 7, ...
      'InsetBox1',[33.28 34.58 -118.5 -118], ...
      'InsetPos1', [0.6 0.25 0.2 0.3], ...
      'InsetBox2',[36.83 38.62 -123.16 -121.84], ...
      'InsetPos2', [0.6 0.6 0.2 0.3]);


