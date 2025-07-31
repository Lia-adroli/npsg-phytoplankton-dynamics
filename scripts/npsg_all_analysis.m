% MASTER SCRIPT FOR NPSG PHYTOPLANKTON AND OHC ANALYSIS
% Author: Lia Adroli (2025)

%% === 1. MIXED LAYER DEPTH ANALYSIS ===
function run_mld_analysis(input_csv)
T = readtable(input_csv);
T.TIME = datetime(T.TIME, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
T.YEAR = year(T.TIME);
T.MONTH = month(T.TIME);
valid = isfinite(T.TEMP_ADJUSTED) & isfinite(T.PSAL_ADJUSTED) & isfinite(T.PRES_ADJUSTED) & T.PRES_ADJUSTED >= 0 & T.PRES_ADJUSTED <= 1000;
T = T(valid, :);
depth_bins = 0:5:1000;
depth_centers = depth_bins(1:end-1) + diff(depth_bins)/2;
all_times = unique(dateshift(datetime(T.YEAR, T.MONTH, 1), 'start', 'month'));
time_list = datetime.empty;
mld_list = [];
for t = 1:length(all_times)
    this_month = all_times(t);
    Ym = T(T.YEAR == year(this_month) & T.MONTH == month(this_month), :);
    if height(Ym) < 10, continue; end
    mean_temp = NaN(size(depth_centers));
    mean_pres = NaN(size(depth_centers));
    for d = 1:length(depth_centers)
        in_bin = Ym.PRES_ADJUSTED >= depth_bins(d) & Ym.PRES_ADJUSTED < depth_bins(d+1);
        if sum(in_bin) >= 5
            mean_temp(d) = mean(Ym.TEMP_ADJUSTED(in_bin), 'omitnan');
            mean_pres(d) = mean(Ym.PRES_ADJUSTED(in_bin), 'omitnan');
        end
    end
    smooth_temp = smoothdata(mean_temp, 'movmean', 3, 'omitnan');
    smooth_pres = smoothdata(mean_pres, 'movmean', 3, 'omitnan');
    if all(isnan(smooth_temp)), continue; end
    ref_idx = find(~isnan(smooth_pres) & smooth_pres >= 5, 1);
    if isempty(ref_idx) || isnan(smooth_temp(ref_idx)), continue; end
    T5 = smooth_temp(ref_idx);
    dT = abs(smooth_temp - T5);
    mld_temp_idx = find(dT > 0.2, 1);
    if ~isempty(mld_temp_idx)
        mld_depth = smooth_pres(mld_temp_idx);
        time_list(end+1, 1) = this_month;
        mld_list(end+1, 1) = mld_depth;
    end
end
MLD_Table = table(time_list, mld_list, 'VariableNames', {'Time', 'MLD_TempOnly'});
writetable(MLD_Table, 'monthly_mld_temperature_only_argo1.csv');
end

%% === 2. TREND ANALYSIS WITH 3x STD OUTLIER REMOVAL ===
function run_trend_analysis(input_csv)
T = readtable(input_csv);
T.time = datetime(T.Date, 'InputFormat', 'dd-MMM-yy');
T.DOY = day(T.time, 'dayofyear');
T.year = year(T.time);
vars = {'SST', 'AOD_500_555_nm', 'Euphotic_Depth', 'SLA', 'Evapotranspiration', 'Cloud_Fraction', 'Surface_Current', 'Precipitation', 'Surface_Wind', 'Picophytoplankton', 'Surface_Chlorophyll', 'Nanophytoplankton', 'Microphytoplankton'};
results = table('Size', [length(vars), 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, 'VariableNames', {'Variable', 'Slope', 'R2', 'P_value'});
for i = 1:length(vars)
    varname = vars{i};
    y = T.(varname);
    good = true(size(y));
    for d = 1:366
        idx = T.DOY == d;
        mu = mean(y(idx), 'omitnan');
        sigma = std(y(idx), 'omitnan');
        good(idx) = good(idx) & abs(y(idx) - mu) <= 3*sigma;
    end
    x = year(T.time(good)) + (day(T.time(good), 'dayofyear') - 1)/365;
    y = y(good);
    valid = isfinite(x) & isfinite(y);
    if sum(valid) >= 30
        mdl = fitlm(x(valid), y(valid));
        results.Variable(i) = string(varname);
        results.Slope(i) = mdl.Coefficients.Estimate(2);
        results.R2(i) = mdl.Rsquared.Adjusted;
        results.P_value(i) = mdl.Coefficients.pValue(2);
    end
end
writetable(results, 'trend_results.csv');
end

%% === 3. CROSS-CORRELATION WITH OUTLIER FILTER ===
function run_crosscorr_analysis(input_csv)
T = readtable(input_csv);
T.time = datetime(T.Date, 'InputFormat', 'dd-MMM-yy');
T.DOY = day(T.time, 'dayofyear');
vars = {'SST', 'AOD_500_555_nm', 'Euphotic_Depth', 'SLA', 'Evapotranspiration', 'Cloud_Fraction', 'Surface_Current', 'Precipitation', 'Surface_Wind', 'Picophytoplankton', 'Surface_Chlorophyll', 'Nanophytoplankton', 'Microphytoplankton'};
for i = 1:length(vars)
    y = T.(vars{i});
    for d = 1:366
        idx = T.DOY == d;
        mu = mean(y(idx), 'omitnan');
        sigma = std(y(idx), 'omitnan');
        y(idx) = y(idx) .* (abs(y(idx)-mu) <= 3*sigma);
    end
    T.(vars{i}) = y;
end
T.MonthID = dateshift(T.time, 'start', 'month');
G = findgroups(T.MonthID);
monthly.Date = splitapply(@(x) x(1), T.MonthID, G);
for i = 1:length(vars)
    monthly.(vars{i}) = splitapply(@(x) mean(x, 'omitnan'), T.(vars{i}), G);
end
ref_var = 'SST';
lags = -12:12;
C = NaN(length(vars), length(lags));
for i = 1:length(vars)
    x = monthly.(ref_var);
    y = monthly.(vars{i});
    valid = isfinite(x) & isfinite(y);
    x = x(valid); y = y(valid);
    x = x - mean(x, 'omitnan');
    y = y - mean(y, 'omitnan');
    [xc, ~] = xcorr(x, y, 12, 'coeff');
    C(i, :) = xc;
end
save('crosscorr_results.mat', 'C', 'vars', 'lags');
end

%% === 4. GRANGER CAUSALITY USING RAW ANOMALIES ===
function run_granger_test(input_csv)
T = readtable(input_csv);
T.time = datetime(T.Date, 'InputFormat', 'dd-MMM-yy');
ref_var = 'Surface_Chlorophyll';
vars = {'SST', 'AOD_500_555_nm', 'Euphotic_Depth', 'SLA', 'Evapotranspiration', 'Cloud_Fraction', 'Surface_Current', 'Precipitation', 'Surface_Wind', 'Picophytoplankton', 'Nanophytoplankton', 'Microphytoplankton'};
max_lag = 30;
for v = [vars, {ref_var}]
    name = v{1};
    data = T.(name);
    T.(name) = (data - mean(data, 'omitnan')) ./ std(data, 'omitnan');
end
results = table('Size', [length(vars), 3], 'VariableTypes', {'string', 'double', 'double'}, 'VariableNames', {'Variable', 'F_stat', 'P_value'});
for i = 1:length(vars)
    x = T.(vars{i});
    y = T.(ref_var);
    valid = isfinite(x) & isfinite(y);
    x = x(valid); y = y(valid);
    if length(x) <= max_lag + 10, continue; end
    Y = y(max_lag+1:end);
    X_lag = lagmatrix([x y], 1:max_lag);
    X_g = X_lag(max_lag+1:end, 1:max_lag);
    X_y = X_lag(max_lag+1:end, max_lag+1:end);
    mask = all(~isnan([Y, X_g, X_y]), 2);
    mdl_full = fitlm([X_y(mask,:) X_g(mask,:)], Y(mask));
    mdl_null = fitlm(X_y(mask,:), Y(mask));
    F = ((mdl_null.SSE - mdl_full.SSE) / max_lag) / (mdl_full.SSE / mdl_full.DFE);
    p = 1 - fcdf(F, max_lag, mdl_full.DFE);
    results.Variable(i) = vars{i};
    results.F_stat(i) = F;
    results.P_value(i) = p;
end
writetable(results, 'granger_results.csv');
end

%% === 5. WAVELET ON RAW DAILY CHLA ===
function run_wavelet_analysis(input_csv)
data = readtable(input_csv, 'VariableNamingRule', 'preserve');
data.Date = datetime(data.Date, 'InputFormat', 'dd/MM/yyyy HH:mm');
chl = data.chl;
valid_idx = ~isnan(chl);
chl = chl(valid_idx);
time = data.Date(valid_idx);
[wt, f] = cwt(chl, 'amor', years(1/365));
f_numeric = years(f);
figure('Position', [100 100 1000 400]);
imagesc(time, f_numeric, abs(wt));
set(gca, 'YDir', 'normal');
xlabel('Time'); ylabel('Period (years)');
title('Wavelet Power Spectrum of Chlorophyll-a');
colormap(flipud(brewermap([], 'RdYlBu'))); colorbar;
end

%% === 6. CROSS-WAVELET CHLA ANOM VS ENSO/PDO ===
function run_crosswavelet_analysis(input_csv)
T = readtable(input_csv);
T.Date = datetime(T.Date);
chl = T.chl_anom;
enso = T.ENSO_index;
pdo = T.PDO_index;
dt = 1/12;
[Wxy1, period1, ~, coi1] = cross_wavelet(chl, enso, dt);
figure;
imagesc(T.Date, log2(period1), abs(Wxy1)); axis xy;
xlabel('Time'); ylabel('Log2 Period'); title('Chla vs ENSO');
hold on; plot(T.Date, log2(coi1), 'w--'); colorbar;
[Wxy2, period2, ~, coi2] = cross_wavelet(chl, pdo, dt);
figure;
imagesc(T.Date, log2(period2), abs(Wxy2)); axis xy;
xlabel('Time'); ylabel('Log2 Period'); title('Chla vs PDO');
hold on; plot(T.Date, log2(coi2), 'w--'); colorbar;
end

%% === 7. SSTa & OHC GRIDDING AND WARMING AREA ===
function run_ssta_ohc_analysis(ssta_nc, ohc_nc)
lat = ncread(ssta_nc, 'latitude');
lon = ncread(ssta_nc, 'longitude');
time = ncread(ssta_nc, 'time');
ssta = ncread(ssta_nc, 'sst_anomaly');
time_dt = datetime(1950,1,1) + days(time);
[nlat, nlon, nt] = size(ssta);
ssta_2d = reshape(ssta, [], nt);
ssta_trend = NaN(size(ssta_2d,1), 1);
for i = 1:size(ssta_2d,1)
    y = ssta_2d(i,:)';
    t = datenum(time_dt - min(time_dt));
    if sum(isfinite(y)) > 24
        mdl = fitlm(t(isfinite(y)), y(isfinite(y)));
        ssta_trend(i) = mdl.Coefficients.Estimate(2) * 365;
    end
end
ssta_trend_map = reshape(ssta_trend, nlat, nlon);
area_map = ones(nlat,nlon);  % Placeholder
expansion_area = NaN(nt,1);
for t = 1:nt
    mask = ssta(:,:,t) > 1;
    expansion_area(t) = nansum(area_map(mask), 'all');
end
figure;
plot(time_dt, expansion_area / 1e6, 'k', 'LineWidth', 2);
xlabel('Year'); ylabel('Area >1°C SSTa (M km²)');
title('Warming Expansion Area');

ohc = ncread(ohc_nc, 'ohc');
ohc_2d = reshape(ohc, [], nt);
ohc_trend = NaN(size(ohc_2d,1), 1);
for i = 1:size(ohc_2d,1)
    y = ohc_2d(i,:)';
    t = datenum(time_dt - min(time_dt));
    if sum(isfinite(y)) > 24
        mdl = fitlm(t(isfinite(y)), y(isfinite(y)));
        ohc_trend(i) = mdl.Coefficients.Estimate(2) * 365;
    end
end
ohc_trend_map = reshape(ohc_trend, nlat, nlon);
end

%% === 8. BGC-ARGO DEPTH-STRATIFIED TRENDS ===
function run_bgc_argo_depth_change(input_csv)
T = readtable(input_csv);
T.time = datetime(T.time);
depth = T.depth;
chl = T.chl;
nitrate = T.nitrate;
oxygen = T.oxygen;
early_idx = T.time < datetime(2021,1,1);
recent_idx = T.time >= datetime(2023,1,1);
zbin = 0:10:200;
zmid = zbin(1:end-1) + 5;
bin_mean = @(x,z) arrayfun(@(k) mean(x(z >= zbin(k) & z < zbin(k+1)), 'omitnan'), 1:length(zmid));
chl_early = bin_mean(chl(early_idx), depth(early_idx));
chl_recent = bin_mean(chl(recent_idx), depth(recent_idx));
n_early = bin_mean(nitrate(early_idx), depth(early_idx));
n_recent = bin_mean(nitrate(recent_idx), depth(recent_idx));
o2_early = bin_mean(oxygen(early_idx), depth(early_idx));
o2_recent = bin_mean(oxygen(recent_idx), depth(recent_idx));
figure;
subplot(1,3,1); plot(chl_early, -zmid, 'b', chl_recent, -zmid, 'r'); xlabel('Chl'); ylabel('Depth'); title('Chl');
subplot(1,3,2); plot(n_early, -zmid, 'b', n_recent, -zmid, 'r'); xlabel('NO₃⁻'); title('Nitrate');
subplot(1,3,3); plot(o2_early, -zmid, 'b', o2_recent, -zmid, 'r'); xlabel('O₂'); title('Oxygen');
sgtitle('Vertical Changes in BGC Profiles');
fprintf('Chl-a DCM change: %.1f%%\n', 100*(max(chl_recent)-max(chl_early))/max(chl_early));
fprintf('NO₃⁻ DCM change: %.1f%%\n', 100*(max(n_recent)-max(n_early))/max(n_early));
fprintf('O₂ DCM change: %.1f%%\n', 100*(max(o2_recent)-max(o2_early))/max(o2_early));
end
