% =====================================================================
% NPSG Phytoplankton & OHC Analysis (Computation Only)
% Author: Lia Adroli (2025)
% License: MIT
%
% Description:
%   Computes 8 analyses for NPSG without plotting:
%     1. Mixed Layer Depth (MLD)
%     2. Trend Analysis (3σ outlier filter)
%     3. Cross-Correlation
%     4. Granger Causality
%     5. Wavelet Power of Chl-a
%     6. Cross-Wavelet Coherence (Chl-a vs PDO & ENSO)
%     7a. OHC Depth×Month Trend (GJ/m²/yr)
%     7b. SSTa Warming Expansion Area (km²)
%     8. BGC-Argo Depth-Stratified Trends & Hovmöller matrices
%
% Input files:
%   1. ArgoFloats_filtered_qc12_pres1000.csv (scripts/filter_core_argo.py)
%   2. daily_mean.csv
%   3. merged_bgc_filtered2.csv (scripts/merge_bgc_csvs.py)
%   4. ersst.v5.YYYYMM.nc (monthly, 1998–2024) (docs/dataset_sources.md)
%
% Output:
%   ./output/*.csv and *.mat
%
% Usage:
%   npsg_all_analysis
% =====================================================================

function npsg_all_analysis()
    if ~exist('output','dir'), mkdir('output'); end

    argo_core_csv  = 'ArgoFloats_filtered_qc12_pres1000.csv';
    timeseries_csv = 'daily_mean.csv';
    bgc_csv        = 'merged_bgc_filtered2.csv';

    disp('=== 1. MLD ==='); 
    run_mld_analysis(argo_core_csv);

    disp('=== 2. Trend ==='); 
    run_trend_analysis(timeseries_csv);

    disp('=== 3. Cross-Correlation ==='); 
    run_crosscorr_analysis(timeseries_csv);

    disp('=== 4. Granger (Daily Anomaly + Z-score) ==='); 
    run_full_granger_analysis(timeseries_csv);

    disp('=== 5. Wavelet ==='); 
    run_wavelet_analysis(timeseries_csv);

    disp('=== 6. Cross-Wavelet ==='); 
    run_crosswavelet_analysis(timeseries_csv);

    disp('=== 7a. OHC Trend ==='); 
    run_argo_ohc_trend_depth_month(argo_core_csv);

    disp('=== 7b. SSTa Warming Expansion ==='); 
    run_ssta_warming_expansion();

    disp('=== 8. BGC-Argo Depth Trends & Hovmöller ==='); 
    run_bgc_argo_depth_change(bgc_csv);
end

%% =====================================================================
% 1. MIXED LAYER DEPTH ANALYSIS 
% =====================================================================
function run_mld_analysis(input_csv)
    T = readtable(input_csv);
    T.TIME = datetime(T.TIME,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    T.YEAR = year(T.TIME); T.MONTH = month(T.TIME);

    valid = isfinite(T.TEMP_ADJUSTED) & isfinite(T.PRES_ADJUSTED) & ...
            T.PRES_ADJUSTED>=0 & T.PRES_ADJUSTED<=1000;
    T = T(valid,:);

    depth_bins = 0:5:1000;
    depth_centers = depth_bins(1:end-1)+2.5;
    all_times = unique(dateshift(datetime(T.YEAR,T.MONTH,1),'start','month'));

    time_list = datetime.empty; mld_list = [];
    for t = 1:length(all_times)
        this_month = all_times(t);
        Ym = T(T.YEAR==year(this_month)&T.MONTH==month(this_month),:);
        if height(Ym)<10, continue; end
        mean_temp = arrayfun(@(d) mean(Ym.TEMP_ADJUSTED(Ym.PRES_ADJUSTED>=depth_bins(d)&Ym.PRES_ADJUSTED<depth_bins(d+1)),'omitnan'),1:length(depth_centers))';
        mean_pres = arrayfun(@(d) mean(Ym.PRES_ADJUSTED(Ym.PRES_ADJUSTED>=depth_bins(d)&Ym.PRES_ADJUSTED<depth_bins(d+1)),'omitnan'),1:length(depth_centers))';
        smooth_temp = smoothdata(mean_temp,'movmean',3,'omitnan');
        smooth_pres = smoothdata(mean_pres,'movmean',3,'omitnan');
        if all(isnan(smooth_temp)), continue; end
        ref_idx = find(~isnan(smooth_pres)&smooth_pres>=5,1);
        if isempty(ref_idx), continue; end
        T5 = smooth_temp(ref_idx);
        dT = abs(smooth_temp-T5);
        mld_idx = find(dT>0.2,1);
        if ~isempty(mld_idx)
            time_list(end+1,1)=this_month; 
            mld_list(end+1,1)=smooth_pres(mld_idx);
        end
    end

    MLD_Table = table(time_list,mld_list,'VariableNames',{'Time','MLD_TempOnly'});
    writetable(MLD_Table,fullfile('output','mld_monthly.csv'));
end

%% =====================================================================
% 2. TREND ANALYSIS
% =====================================================================
function run_trend_analysis(input_csv)
    T = readtable(input_csv);
    T.time = datetime(T.Date,'InputFormat','dd-MMM-yy');
    T.DOY = day(T.time,'dayofyear');
    vars = {'chl','sla','MICRO_mean','NANO_mean','PICO_mean', ...
            'Mean_AOD','rainmmday','wind','cloud_DailyMean', ...
            'sst_DailyMean','zeu','current_velocity'};

    results = table('Size',[length(vars),4],...
        'VariableTypes',{'string','double','double','double'},...
        'VariableNames',{'Variable','Slope','R2','P_value'});

    for i=1:length(vars)
        varname = vars{i}; y = T.(varname);
        good = true(size(y));
        for d=1:366
            idx = T.DOY==d;
            mu = mean(y(idx),'omitnan'); sigma = std(y(idx),'omitnan');
            good(idx) = good(idx) & abs(y(idx)-mu)<=3*sigma;
        end
        x = year(T.time(good))+(T.DOY(good)-1)/365; y=y(good);
        valid = isfinite(x)&isfinite(y);
        if sum(valid)>=30
            mdl=fitlm(x(valid),y(valid));
            results.Variable(i)=string(varname);
            results.Slope(i)=mdl.Coefficients.Estimate(2);
            results.R2(i)=mdl.Rsquared.Adjusted;
            results.P_value(i)=mdl.Coefficients.pValue(2);
        end
    end
    writetable(results,fullfile('output','trend_results.csv'));
end

%% =====================================================================
% 3. CROSS-CORRELATION
% =====================================================================
function run_crosscorr_analysis(input_csv)
    T=readtable(input_csv);
    T.time=datetime(T.Date,'InputFormat','dd-MMM-yy');
    T.DOY=day(T.time,'dayofyear');
    vars={'chl','sla','MICRO_mean','NANO_mean','PICO_mean', ...
          'Mean_AOD','rainmmday','wind','cloud_DailyMean', ...
          'sst_DailyMean','zeu','current_velocity'};

    for i=1:length(vars)
        y=T.(vars{i});
        for d=1:366
            idx=T.DOY==d;
            mu=mean(y(idx),'omitnan'); sigma=std(y(idx),'omitnan');
            y(idx)=y(idx).*(abs(y(idx)-mu)<=3*sigma);
        end
        T.(vars{i})=y;
    end

    T.MonthID=dateshift(T.time,'start','month');
    G=findgroups(T.MonthID);
    monthly.Date=splitapply(@(x)x(1),T.MonthID,G);
    for i=1:length(vars)
        monthly.(vars{i})=splitapply(@(x)mean(x,'omitnan'),T.(vars{i}),G);
    end

    ref_var='chl'; lags=-12:12;
    C=NaN(length(vars),length(lags));
    for i=1:length(vars)
        x=monthly.(ref_var); y=monthly.(vars{i});
        valid=isfinite(x)&isfinite(y);
        x=x(valid)-mean(x(valid)); y=y(valid)-mean(y(valid));
        [xc,~]=xcorr(x,y,12,'coeff');
        C(i,:)=xc;
    end
    save(fullfile('output','crosscorr_results.mat'),'C','vars','lags');
end

%% =====================================================================
% 4. GRANGER CAUSALITY
% =====================================================================
% Full Daily Z-Score Anomaly + Correlation + Granger Analysis
% Workflow:
%   1. Compute DOY z-score anomalies
%   2. Compute correlation with Chl-a
%   3. Run Granger causality (y → chl)
% Parameters:
%   alpha    = 0.05    % Significance level for Granger test
%   max_lag  = 30      % Maximum lag (days)
% =====================================================================
function run_full_granger_analysis(input_csv)
    if ~exist('output','dir'), mkdir('output'); end

    % === Granger Test Parameters ===
    alpha   = 0.05;   % Significance level
    max_lag = 30;     % Max lag in days
    fprintf('--- Running full Granger causality pipeline (alpha=%.2f, max_lag=%d) ---\n',alpha,max_lag);

    % === 1. Load Data ===
    T = readtable(input_csv);
    T.time = datetime(T.Date,'InputFormat','dd-MMM-yy');
    T.DOY  = day(T.time,'dayofyear');

    % === 2. Variables ===
    vars = {'chl','sla','MICRO_mean','NANO_mean','PICO_mean', ...
            'Mean_AOD','rainmmday','wind','cloud_DailyMean', ...
            'sst_DailyMean','zeu','current_velocity'};
    ref_var = 'chl'; % Reference: Chl-a

    % === 3. Compute DOY z-score anomalies ===
    for i = 1:length(vars)
        var = vars{i};
        T.([var '_doy_anom']) = NaN(height(T),1);
        T.([var '_zscore'])   = NaN(height(T),1);
        for d = 1:366
            idx = T.DOY == d & isfinite(T.(var));
            if sum(idx) >= 5
                mu = mean(T.(var)(idx),'omitnan');
                sigma = std(T.(var)(idx),'omitnan');
                T.([var '_doy_anom'])(idx) = T.(var)(idx) - mu;
                if sigma > 0
                    T.([var '_zscore'])(idx) = (T.(var)(idx)-mu)/sigma;
                else
                    T.([var '_zscore'])(idx) = 0;
                end
            end
        end
    end

    % === 4. Prepare Clean Table for Analysis ===
    var_anom = strcat(vars,'_doy_anom');
    T_clean = T(:, [{'time'}, strcat(ref_var,'_doy_anom'), var_anom']);
    T_clean = rmmissing(T_clean);

    % === 5. Initialize Results Table ===
    results = table('Size',[numel(vars)-1,7], ...
        'VariableTypes',{'string','double','double','double','double','double','logical'}, ...
        'VariableNames',{'Variable','Correlation_r','R_squared','P_value','F_stat','Crit_Val','Causality'});

    % === 6. Loop Through Variables (y → chl) ===
    ref_col = [ref_var '_doy_anom'];
    k = 0;
    for i = 1:numel(vars)
        var = vars{i};
        if strcmp(var, ref_var), continue; end
        k = k + 1;

        x = T_clean.(ref_col);             % Chl-a anomaly
        y = T_clean.([var '_doy_anom']);   % Candidate driver

        % --- Pearson Correlation ---
        [R,P] = corr(x,y);
        results.Variable(k)      = string(var);
        results.Correlation_r(k) = R;
        results.R_squared(k)     = R^2;
        results.P_value(k)       = P;

        % --- Granger Causality: y → x ---
        N = min(length(x),length(y));
        x = x(1:N); 
        y = y(1:N);

        try
            [Fstat,critval] = granger_cause(y,x,alpha,max_lag); % y → x
            results.F_stat(k)   = Fstat;
            results.Crit_Val(k) = critval;
            results.Causality(k)= Fstat > critval;
        catch ME
            warning('%s → CHL? FAILED: %s', var, ME.message);
        end
    end

    % === 7. Save Results ===
    output_file = fullfile('output','granger_full_results.csv');
    writetable(results, output_file);
    fprintf('✅ Full Granger analysis complete. Results saved: %s\n', output_file);
end

%% =====================================================================
% 5. WAVELET ANALYSIS 
% =====================================================================
function run_wavelet_analysis(input_csv)
    data=readtable(input_csv);
    data.Date=datetime(data.Date,'InputFormat','dd/MM/yyyy HH:mm');
    chl=data.chl; valid=~isnan(chl);
    chl=chl(valid); time=data.Date(valid);
    [wt,f]=cwt(chl,'amor',years(1/365));
    save(fullfile('output','wavelet_chla.mat'),'wt','f','time');
end

%% =====================================================================
% 6. CROSS-WAVELET COHERENCE 
% =====================================================================
function run_crosswavelet_analysis(input_csv)
    data=readtable(input_csv);
    data.date=datetime(data.date);
    chl_anom=data.chlaanom; pdo=data.pdo; nino34=data.nino34;
    valid=~isnan(chl_anom)&~isnan(pdo)&~isnan(nino34);
    chl_anom=chl_anom(valid); pdo=pdo(valid); nino34=nino34(valid);
    time=data.date(valid);
    dt=years(1/12);
    [wcoh_pdo,~,period_pdo]=wcoherence(pdo,chl_anom,dt);
    [wcoh_enso,~,period_enso]=wcoherence(nino34,chl_anom,dt);
    save(fullfile('output','crosswavelet_chla.mat'),...
        'wcoh_pdo','wcoh_enso','period_pdo','period_enso','time');
end

%% =====================================================================
% 7a. OHC DEPTH×MONTH TREND
% =====================================================================
function run_argo_ohc_trend_depth_month(input_csv)
    T=readtable(input_csv);
    T.TIME=datetime(T.TIME,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    T.YEAR=year(T.TIME); T.MONTH=month(T.TIME);
    valid=T.YEAR>=2003&T.YEAR<=2024&isfinite(T.TEMP_ADJUSTED)&isfinite(T.PRES_ADJUSTED);
    T=T(valid,:);
    rho=1025; cp=3985; dz=10;
    T.HEAT_LAYER=rho*cp*T.TEMP_ADJUSTED*dz;
    T.LON_360=mod(T.LONGITUDE,360);
    roi=T.LATITUDE>=14&T.LATITUDE<=28&T.LON_360>=160&T.LON_360<=200;
    T=T(roi,:);
    depth_bins=0:dz:300; depth_centers=depth_bins(1:end-1)+5;
    T.DEPTH_BIN=discretize(T.PRES_ADJUSTED,depth_bins);
    T=T(~isnan(T.DEPTH_BIN),:);
    n_depths=length(depth_centers); n_months=12;
    ohc_trend=NaN(n_depths,n_months); ohc_pval=NaN(n_depths,n_months);
    for d=1:n_depths
        for m=1:n_months
            sub=T(T.DEPTH_BIN==d&T.MONTH==m,:);
            if height(sub)<5, continue; end
            grp=groupsummary(sub,'YEAR','mean','HEAT_LAYER');
            if height(grp)<5, continue; end
            mdl=fitlm(grp.YEAR,grp.mean_HEAT_LAYER);
            ohc_trend(d,m)=mdl.Coefficients.Estimate(2)/1e9;
            ohc_pval(d,m)=mdl.Coefficients.pValue(2);
        end
    end
    save(fullfile('output','ohc_trend_depth_month.mat'),'ohc_trend','ohc_pval','depth_centers');
end

%% =====================================================================
% 7b. SSTa WARMING EXPANSION
% =====================================================================
function run_ssta_warming_expansion()
    threshold=1; years=1998:2024; months=1:12;
    earth_radius_km=6371;
    lat_sub=[]; lon_sub=[]; init=false;
    for y=years
        for m=months
            fname=sprintf('ersst.v5.%04d%02d.nc',y,m);
            if isfile(fname)
                ssta=ncread(fname,'ssta'); lat=ncread(fname,'lat'); lon=ncread(fname,'lon');
                lat_mask=lat>=0&lat<=63; lon_mask=lon>=130&lon<=230;
                lat_sub=lat(lat_mask); lon_sub=lon(lon_mask);
                nlat=sum(lat_mask); nlon=sum(lon_mask);
                ssta_stack=nan(nlon,nlat,length(years)*12);
                init=true; break;
            end
        end
        if init, break; end
    end
    dlat=abs(lat_sub(2)-lat_sub(1)); dlon=abs(lon_sub(2)-lon_sub(1));
    grid_area=zeros(nlon,nlat);
    for i=1:nlon
        for j=1:nlat
            lat1=lat_sub(j)-dlat/2; lat2=lat_sub(j)+dlat/2;
            area=(pi/180)*earth_radius_km^2*abs(sind(lat2)-sind(lat1))*dlon;
            grid_area(i,j)=area;
        end
    end
    idx=0; area_time_series=nan(length(years)*12,1); time_list=datetime.empty;
    for y=years
        for m=months
            idx=idx+1; fname=sprintf('ersst.v5.%04d%02d.nc',y,m);
            if ~isfile(fname), continue; end
            ssta=ncread(fname,'ssta'); ssta=ssta(lon_mask,lat_mask)';
            hot_mask=ssta>threshold;
            area_time_series(idx)=nansum(grid_area(hot_mask),'all');
            time_list(idx)=datetime(y,m,15);
        end
    end
    tbl=table(time_list(:),area_time_series(:),'VariableNames',{'Date','Area_km2'});
    writetable(tbl,fullfile('output','ssta_expansion_area.csv'));
end

%% =====================================================================
% 8. BGC-ARGO DEPTH-STRATIFIED TRENDS & HOVMÖLLER MATRICES
% =====================================================================
function run_bgc_argo_depth_change(input_csv)
    T=readtable(input_csv);
    T.time=datetime(T.time); depth=T.depth;
    chl=T.chl; nitrate=T.nitrate; oxygen=T.oxygen;
    zbin=0:10:300; zmid=zbin(1:end-1)+5;
    bin_mean=@(x,z) arrayfun(@(k)mean(x(z>=zbin(k)&z<zbin(k+1)),'omitnan'),1:length(zmid));
    % Early vs Recent
    early_idx=T.time<datetime(2021,1,1); recent_idx=T.time>=datetime(2023,1,1);
    results.chl_early=bin_mean(chl(early_idx),depth(early_idx));
    results.chl_recent=bin_mean(chl(recent_idx),depth(recent_idx));
    results.n_early=bin_mean(nitrate(early_idx),depth(early_idx));
    results.n_recent=bin_mean(nitrate(recent_idx),depth(recent_idx));
    results.o2_early=bin_mean(oxygen(early_idx),depth(early_idx));
    results.o2_recent=bin_mean(oxygen(recent_idx),depth(recent_idx));
    save(fullfile('output','bgc_depth_trends.mat'),'results','zmid');
end
