function results = npsg_cbpm_export_pipeline(datadir, opts)
% NPSG_CBPM_EXPORT_PIPELINE
% One-file GitHub-ready function that:
%  - Reads Argo core & BGC, monthly Kd/Zeu/PAR, MLD
%  - Bins to 10 m, computes CbPM NPP, AOU -> AOUc
%  - Computes isopycnal correlations & export proxies
%  - Returns Particle Export Efficiency (PEE, %) for 2020 & 2024
%
% USAGE:
%   results = npsg_cbpm_export_pipeline('../data', struct('save_csv', true));
%
% INPUTS:
%   datadir  : folder with CSVs
%   opts     : (optional) struct with fields:
%       .core_file   (default 'ArgoFloats_filtered_qc12_pres1000.csv')
%       .bgc_file    (default 'merged_bgc_filtered5.csv')
%       .file_kd     (default 'monthly_Kd490_Zeu_ZSD.csv')
%       .file_mld    (default 'monthly_mld_temperature_only_argo1_2003.csv')
%       .years       (default [2020 2024])
%       .months      (default 1:11)  % Jan–Nov
%       .save_csv    (default false) % writes results/pee_summary.csv if true
%
% OUTPUT:
%   results : struct with fields:
%       .PEE_2020_percent, .PEE_2024_percent
%       .PEE_table (2x2 table: Year, PEE_percent)
%       .lag_best_2020, .lag_best_2024
%       .T20_iso, .T24_iso    (isopycnal correlation summaries)
%       .arrays               (NPP/AOUc profiles, masks, KPAR, etc.)
%
% DEPENDENCIES (optional but supported):
%   GSW or SEAWATER toolboxes (automatic fallback included if missing)

% --------------------------- defaults ---------------------------
if nargin < 1 || isempty(datadir), datadir = './data'; end
if nargin < 2, opts = struct; end
getOpt = @(s,f,def) (isfield(s,f) && ~isempty(s.(f))) .* s.(f) + ...
                    (~(isfield(s,f) && ~isempty(s.(f)))) .* def;

core_file = fullfile(datadir, getOpt(opts,'core_file','ArgoFloats_filtered_qc12_pres1000.csv'));
bgc_file  = fullfile(datadir, getOpt(opts,'bgc_file','merged_bgc_filtered5.csv'));
file_kd   = fullfile(datadir, getOpt(opts,'file_kd','monthly_Kd490_Zeu_ZSD.csv'));
file_mld  = fullfile(datadir, getOpt(opts,'file_mld','monthly_mld_temperature_only_argo1_2003.csv'));
YEARS     = getOpt(opts,'years',[2020 2024]);
MONTHS    = getOpt(opts,'months',1:11);
save_csv  = logical(getOpt(opts,'save_csv',false));

repo  = fileparts(fileparts(mfilename('fullpath')));
resdir = fullfile(repo,'results'); if ~exist(resdir,'dir'), mkdir(resdir); end

% constants
depth_bins    = 0:10:300;           
depth_centers = depth_bins(1:end-1) + diff(depth_bins)/2;
rho = 1025;                         % kg m^-3
O2_to_C = 117/170;                  % Redfield C:O2
gC_to_molC = 1/12.011;              % mol C per g C
dt_days = 30.4;

% ---------------------- 1) CORE (T/S/P) ------------------------
T_core = readtable(core_file);
T_core.TIME  = datetime(T_core.TIME, 'InputFormat','yyyy-MM-dd HH:mm:ss');
T_core.YEAR  = year(T_core.TIME);
T_core.MONTH = month(T_core.TIME);
valid_core = isfinite(T_core.TEMP_ADJUSTED) & isfinite(T_core.PSAL_ADJUSTED) & ...
             isfinite(T_core.PRES_ADJUSTED) & T_core.PRES_ADJUSTED <= 300;
T_core = T_core(valid_core, :);

% ---------------------- 2) BGC (DOXY/CHL/BBP) ------------------
T_bgc = readtable(bgc_file);
if iscell(T_bgc.TIME)
    T_bgc.TIME = datetime(T_bgc.TIME,'InputFormat','dd-MMM-yyyy');
else
    T_bgc.TIME = datetime(T_bgc.TIME);
end
T_bgc.YEAR  = year(T_bgc.TIME);
T_bgc.MONTH = month(T_bgc.TIME);
valid_bgc = isfinite(T_bgc.PRES) & T_bgc.PRES <= 300;
T_bgc = T_bgc(valid_bgc, :);

% ---------------------- 3) LIGHT CLEANING -----------------------
if ismember('NITRATE_ADJUSTED', T_bgc.Properties.VariableNames)
    T_bgc.NITRATE_ADJUSTED(T_bgc.NITRATE_ADJUSTED > 20) = NaN;
end
if ismember('CHLA_ADJUSTED', T_bgc.Properties.VariableNames)
    T_bgc.CHLA_ADJUSTED(T_bgc.CHLA_ADJUSTED < 0) = NaN;
end
if ismember('DOXY_ADJUSTED', T_bgc.Properties.VariableNames)
    T_bgc.DOXY_ADJUSTED(T_bgc.DOXY_ADJUSTED < 0 | T_bgc.DOXY_ADJUSTED > 300) = NaN;
end
if ismember('BBP700_ADJUSTED', T_bgc.Properties.VariableNames)
    T_bgc.BBP700_ADJUSTED(T_bgc.BBP700_ADJUSTED < 0 | T_bgc.BBP700_ADJUSTED > 0.05) = NaN;
end

% ---------------------- 4) BIN 10 m, MONTHLY --------------------
T_core.DEPTH_BIN = discretize(T_core.PRES_ADJUSTED, depth_bins);
T_core = T_core(~isnan(T_core.DEPTH_BIN), :);
core_keep = intersect({'TEMP_ADJUSTED','PSAL_ADJUSTED','PRES_ADJUSTED'}, ...
                      T_core.Properties.VariableNames, 'stable');
core_use  = unique([{'YEAR','MONTH','DEPTH_BIN'}, core_keep], 'stable');
T_core    = T_core(:, core_use);
core_bin  = groupsummary(T_core, {'YEAR','MONTH','DEPTH_BIN'}, 'mean', core_keep);
for k = 1:numel(core_keep)
    old = "mean_" + core_keep{k};
    if ismember(old, core_bin.Properties.VariableNames)
        core_bin.Properties.VariableNames{old} = core_keep{k};
    end
end
core_bin.DEPTH_M = depth_centers(core_bin.DEPTH_BIN)';

T_bgc.DEPTH_BIN = discretize(T_bgc.PRES, depth_bins);
T_bgc = T_bgc(~isnan(T_bgc.DEPTH_BIN), :);
bgc_keep = intersect({'CHLA_ADJUSTED','NITRATE_ADJUSTED','BBP700_ADJUSTED','DOXY_ADJUSTED'}, ...
                     T_bgc.Properties.VariableNames, 'stable');
bgc_use  = unique([{'YEAR','MONTH','DEPTH_BIN'}, bgc_keep], 'stable');
T_bgc    = T_bgc(:, bgc_use);
if isempty(bgc_keep)
    bgc_bin = unique(T_bgc(:, {'YEAR','MONTH','DEPTH_BIN'}), 'rows');
else
    bgc_bin = groupsummary(T_bgc, {'YEAR','MONTH','DEPTH_BIN'}, 'mean', bgc_keep);
    for k = 1:numel(bgc_keep)
        old = "mean_" + bgc_keep{k};
        if ismember(old, bgc_bin.Properties.VariableNames)
            bgc_bin.Properties.VariableNames{old} = bgc_keep{k};
        end
    end
end
bgc_bin.DEPTH_M = depth_centers(bgc_bin.DEPTH_BIN)';

% ---------------------- 5) MERGE CORE + BGC --------------------
merge_keys = {'YEAR','MONTH','DEPTH_BIN'};
all_bin = outerjoin(core_bin, bgc_bin, 'Keys', merge_keys, 'MergeKeys', true, 'Type','full');
if ~ismember('DEPTH_M', all_bin.Properties.VariableNames)
    all_bin.DEPTH_M = depth_centers(all_bin.DEPTH_BIN)';
end
all_bin = sortrows(all_bin, {'YEAR','MONTH','DEPTH_BIN'});

% ---------------------- 6) AOU (Garcia & Gordon 1992) ----------
hasT  = ismember('TEMP_ADJUSTED', all_bin.Properties.VariableNames);
hasS  = ismember('PSAL_ADJUSTED', all_bin.Properties.VariableNames);
hasDO = ismember('DOXY_ADJUSTED', all_bin.Properties.VariableNames);
if hasT && hasS && hasDO
    tC = all_bin.TEMP_ADJUSTED;  S = all_bin.PSAL_ADJUSTED;  DO = all_bin.DOXY_ADJUSTED; % µmol/kg
    Ts = log((298.15 - tC) ./ (273.15 + tC));  % scaled temp
    A0=5.80871; A1=3.20291; A2=4.17887; A3=5.10006; A4=-9.86643e-2; A5=3.80369;
    B0=-7.01577e-3; B1=-7.70028e-3; B2=-1.13864e-2; B3=-9.51519e-3; C0=-2.75915e-7;
    lnC = A0 + A1.*Ts + A2.*Ts.^2 + A3.*Ts.^3 + A4.*Ts.^4 + A5.*Ts.^5 + ...
          S.*(B0 + B1.*Ts + B2.*Ts.^2 + B3.*Ts.^3) + C0.*S.^2;
    Cstar = exp(lnC);                           % µmol/kg
    all_bin.O2_SOL_UMOLKG = Cstar;
    all_bin.AOU_UMOLKG    = Cstar - DO;        % µmol/kg
end

% ---------------------- 7) FILTER YEARS/MONTHS -----------------
all_bin = all_bin(ismember(all_bin.YEAR,YEARS) & ismember(all_bin.MONTH,MONTHS), :);

% ---------------------- 8) Kd/PAR (monthly) --------------------
opts_kd = detectImportOptions(file_kd);
opts_kd = setvartype(opts_kd,'Date','datetime');
opts_kd = setvaropts(opts_kd,'Date','InputFormat','d/MM/uuuu');
numCols = intersect({'Kd490','Zeu_from_Kd490','ZSD','PAR_surf'}, opts_kd.VariableNames);
if ~isempty(numCols), opts_kd = setvartype(opts_kd, numCols, 'double'); end
Tkd = readtable(file_kd, opts_kd);
if ~isdatetime(Tkd.Date), Tkd.Date = datetime(Tkd.Date,'InputFormat','d/MM/uuuu'); end
Tkd.Year  = year(Tkd.Date);  Tkd.Month = month(Tkd.Date);
Gkd = groupsummary(Tkd, {'Year','Month'}, 'median', {'Kd490','Zeu_from_Kd490','PAR_surf'});
get11 = @(Y,v) arrayfun(@(m) median(Gkd.(['median_' v])(Gkd.Year==Y & Gkd.Month==m), 'omitnan'), MONTHS)';

Kd20  = get11(2020,'Kd490');           Kd24  = get11(2024,'Kd490');
ZEU20 = get11(2020,'Zeu_from_Kd490');  ZEU24 = get11(2024,'Zeu_from_Kd490');
PAR20 = get11(2020,'PAR_surf');        PAR24 = get11(2024,'PAR_surf');

% ---------------------- 9) MLD (monthly) -----------------------
opts_mld = detectImportOptions(file_mld);
opts_mld = setvartype(opts_mld, 'Time', 'char');
Tmld = readtable(file_mld, opts_mld);
Tmld.Time = datetime(strtrim(Tmld.Time), 'InputFormat','d/M/uuuu');
bad = isnat(Tmld.Time);
if any(bad), Tmld.Time(bad) = datetime(strtrim(Tmld.Time(bad)), 'InputFormat','dd/MM/uuuu'); end
assert(ismember('MLD_TempOnly', Tmld.Properties.VariableNames), 'Need MLD_TempOnly');
Tmld.MLD_TempOnly = double(Tmld.MLD_TempOnly);
Tmld.Year  = year(Tmld.Time); Tmld.Month = month(Tmld.Time);
MLD20 = arrayfun(@(m) mean(Tmld.MLD_TempOnly(Tmld.Year==2020 & Tmld.Month==m), 'omitnan'), MONTHS)';
MLD24 = arrayfun(@(m) mean(Tmld.MLD_TempOnly(Tmld.Year==2024 & Tmld.Month==m), 'omitnan'), MONTHS)';
MLD20(~isfinite(MLD20)) = nanmedian(MLD20);
MLD24(~isfinite(MLD24)) = nanmedian(MLD24);

% ---------------------- 10) Build grids ------------------------
if ismember('DEPTH_M', all_bin.Properties.VariableNames)
    zc = unique(all_bin.DEPTH_M); zc = zc(~isnan(zc)); zc = sort(zc);
else
    zc = depth_centers(:);
end
nZ = numel(zc); dz = median(diff(zc));
BBP7 = nan(nZ, numel(MONTHS), 2);
CHL  = nan(nZ, numel(MONTHS), 2);

for yi = 1:2
    yy = YEARS(yi);
    dyy = all_bin(all_bin.YEAR==yy & ismember(all_bin.MONTH,MONTHS), :);
    if ~ismember('DEPTH_M', dyy.Properties.VariableNames) && ismember('DEPTH_BIN', dyy.Properties.VariableNames)
        dyy.DEPTH_M = depth_centers(dyy.DEPTH_BIN);
    end
    for mi = 1:numel(MONTHS)
        m = MONTHS(mi); dm = dyy(dyy.MONTH==m, :);
        if isempty(dm), continue; end
        bbp_here = nan(nZ,1); chl_here = nan(nZ,1);
        if ismember('BBP700_ADJUSTED', dm.Properties.VariableNames) && numel(unique(dm.DEPTH_M)) > 1
            bbp_here = interp1(dm.DEPTH_M, dm.BBP700_ADJUSTED, zc, 'linear', NaN);
        elseif ismember('BBP700_ADJUSTED', dm.Properties.VariableNames) && ~isempty(dm.DEPTH_M)
            [~,iz] = min(abs(zc - dm.DEPTH_M(1))); bbp_here(iz) = dm.BBP700_ADJUSTED(1);
        end
        if ismember('CHLA_ADJUSTED', dm.Properties.VariableNames) && numel(unique(dm.DEPTH_M)) > 1
            chl_here = interp1(dm.DEPTH_M, dm.CHLA_ADJUSTED, zc, 'linear', NaN);
        elseif ismember('CHLA_ADJUSTED', dm.Properties.VariableNames) && ~isempty(dm.DEPTH_M)
            [~,iz] = min(abs(zc - dm.DEPTH_M(1))); chl_here(iz) = dm.CHLA_ADJUSTED(1);
        end
        BBP7(:, mi, yi) = bbp_here;
        CHL(:,  mi, yi) = chl_here;
    end
end

% ---------------------- 11) KPAR from Kd, MLD ------------------
KPAR20 = nan(numel(MONTHS),1); KPAR24 = nan(numel(MONTHS),1);
for mi = 1:numel(MONTHS)
    if isfinite(Kd20(mi)) && Kd20(mi) > 0
        invKd = 1/Kd20(mi);
        KPAR20(mi) = (MLD20(mi) <= invKd) * (0.0864 + 0.884*Kd20(mi) - 0.00137*invKd) + ...
                     (MLD20(mi)  > invKd) * (0.0665 + 0.874*Kd20(mi) - 0.00121*invKd);
    end
    if isfinite(Kd24(mi)) && Kd24(mi) > 0
        invKd = 1/Kd24(mi);
        KPAR24(mi) = (MLD24(mi) <= invKd) * (0.0864 + 0.884*Kd24(mi) - 0.00137*invKd) + ...
                     (MLD24(mi)  > invKd) * (0.0665 + 0.874*Kd24(mi) - 0.00121*invKd);
    end
end

% ---------------------- 12) NPP(z,m): g C m^-3 d^-1 -------------
NPP_gm3_d_20 = nan(nZ, numel(MONTHS));
NPP_gm3_d_24 = nan(nZ, numel(MONTHS));
epsC = 1e-6;
for mi = 1:numel(MONTHS)
    Iz20 = PAR20(mi) .* exp(-KPAR20(mi) .* zc); Iz20 = max(Iz20(:),0);
    Iz24 = PAR24(mi) .* exp(-KPAR24(mi) .* zc); Iz24 = max(Iz24(:),0);
    % 2020
    bbp470  = BBP7(:,mi,1) .* (470/700)^(-0.78);
    Cphy_mg = 12128 .* bbp470 + 0.59;
    chl_mg  = CHL(:,mi,1);
    good    = isfinite(Cphy_mg) & isfinite(chl_mg) & (Cphy_mg > epsC);
    ratio   = nan(nZ,1); ratio(good) = chl_mg(good)./Cphy_mg(good);
    mu20    = nan(nZ,1);
    mu20(good) = (2 .* ratio(good) .* (1 - exp(-5 .* Iz20(good)))) ./ ...
                 (0.022 + (0.045-0.022) .* exp(-3 .* Iz20(good)));
    NPP_gm3_d_20(:,mi) = (Cphy_mg .* mu20) ./ 1000;
    % 2024
    bbp470  = BBP7(:,mi,2) .* (470/700)^(-0.78);
    Cphy_mg = 12128 .* bbp470 + 0.59;
    chl_mg  = CHL(:,mi,2);
    good    = isfinite(Cphy_mg) & isfinite(chl_mg) & (Cphy_mg > epsC);
    ratio   = nan(nZ,1); ratio(good) = chl_mg(good)./Cphy_mg(good);
    mu24    = nan(nZ,1);
    mu24(good) = (2 .* ratio(good) .* (1 - exp(-5 .* Iz24(good)))) ./ ...
                 (0.022 + (0.045-0.022) .* exp(-3 .* Iz24(good)));
    NPP_gm3_d_24(:,mi) = (Cphy_mg .* mu24) ./ 1000;
end

% ---------------------- 13) AOU(z,m) arrays ---------------------
AOU20 = nan(nZ, numel(MONTHS));
AOU24 = nan(nZ, numel(MONTHS));
for mi = 1:numel(MONTHS)
    m = MONTHS(mi);
    d20 = all_bin(all_bin.YEAR==2020 & all_bin.MONTH==m, :);
    d24 = all_bin(all_bin.YEAR==2024 & all_bin.MONTH==m, :);
    if ~isempty(d20)
        if ~ismember('DEPTH_M', d20.Properties.VariableNames) && ismember('DEPTH_BIN', d20.Properties.VariableNames)
            d20.DEPTH_M = depth_centers(d20.DEPTH_BIN);
        end
        AOU20(:,mi) = interp1(d20.DEPTH_M, d20.AOU_UMOLKG, zc, 'linear', NaN);
    end
    if ~isempty(d24)
        if ~ismember('DEPTH_M', d24.Properties.VariableNames) && ismember('DEPTH_BIN', d24.Properties.VariableNames)
            d24.DEPTH_M = depth_centers(d24.DEPTH_BIN);
        end
        AOU24(:,mi) = interp1(d24.DEPTH_M, d24.AOU_UMOLKG, zc, 'linear', NaN);
    end
end

% ---------------------- 14) ML AOU & Areal NPP -----------------
ML_AOU20 = nan(numel(MONTHS),1);
ML_AOU24 = nan(numel(MONTHS),1);
for mi = 1:numel(MONTHS)
    ML_AOU20(mi) = mean(AOU20(zc <= MLD20(mi), mi), 'omitnan');
    ML_AOU24(mi) = mean(AOU24(zc <= MLD24(mi), mi), 'omitnan');
end
NPP_areal_20 = nansum(NPP_gm3_d_20 .* dz, 1)';    % g C m^-2 d^-1
NPP_areal_24 = nansum(NPP_gm3_d_24 .* dz, 1)';

% ---------------------- 15) Carbon-unit conversions -------------
AOU_O2_molm3_20 = (AOU20 .* rho) ./ 1e6;                   
AOU_O2_molm3_24 = (AOU24 .* rho) ./ 1e6;
AOU_C_molm3_20  = O2_to_C .* AOU_O2_molm3_20;              
AOU_C_molm3_24  = O2_to_C .* AOU_O2_molm3_24;
NPP_molC_m3_d_20 = NPP_gm3_d_20 .* gC_to_molC;             
NPP_molC_m3_d_24 = NPP_gm3_d_24 .* gC_to_molC;
NPP_molC_m2_d_20 = nansum(NPP_molC_m3_d_20 .* dz, 1)';     
NPP_molC_m2_d_24 = nansum(NPP_molC_m3_d_24 .* dz, 1)';

ZEU20_med = nanmedian(ZEU20); ZEU24_med = nanmedian(ZEU24);
mask_eu20 = zc <= ZEU20_med; mask_eu24 = zc <= ZEU24_med;
AOU_C_belowEZ_m2_20 = nansum(AOU_C_molm3_20(~mask_eu20,:) .* dz, 1)';   
AOU_C_belowEZ_m2_24 = nansum(AOU_C_molm3_24(~mask_eu24,:) .* dz, 1)';

% ---------------------- 16) sigma0 (EOS-80 fallback) -----------
sigma0_eos80 = @(S,T) ( ...
    (999.842594 + 6.793952e-2*T - 9.09529e-3*T.^2 + 1.001685e-4*T.^3 ...
     - 1.120083e-6*T.^4 + 6.536332e-9*T.^5) + ...
    (8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T.^2 - 8.2467e-7*T.^3 + 5.3875e-9*T.^4).*S + ...
    (-5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T.^2).*S.^1.5 + 4.8314e-4*S.^2 ) - 1000;

Tz20 = nan(nZ, numel(MONTHS)); Sz20 = Tz20; sig20 = Tz20;
Tz24 = nan(nZ, numel(MONTHS)); Sz24 = Tz24; sig24 = Tz24;

for yi = 1:2
    YY = YEARS(yi);
    for mi = 1:numel(MONTHS)
        m  = MONTHS(mi);
        tt = all_bin(all_bin.YEAR==YY & all_bin.MONTH==m, :);
        if isempty(tt), continue; end
        if ~ismember('DEPTH_M', tt.Properties.VariableNames) && ismember('DEPTH_BIN', tt.Properties.VariableNames)
            tt.DEPTH_M = depth_centers(tt.DEPTH_BIN);
        end
        Tz = interp1(tt.DEPTH_M, tt.TEMP_ADJUSTED, zc, 'linear', NaN);
        Sz = interp1(tt.DEPTH_M, tt.PSAL_ADJUSTED, zc, 'linear', NaN);
        sg = sigma0_eos80(Sz, Tz);
        if yi==1, Tz20(:,mi)=Tz; Sz20(:,mi)=Sz; sig20(:,mi)=sg;
        else,     Tz24(:,mi)=Tz; Sz24(:,mi)=Sz; sig24(:,mi)=sg; end
    end
end

% ---------------------- 17) Isopycnal correlations -------------
sig_all = [sig20(:); sig24(:)];
sig_all = sig_all(isfinite(sig_all));
sig_min = floor(min(sig_all)*10)/10;
sig_max = ceil( max(sig_all)*10)/10;
sig_edges = sig_min:0.10:sig_max;
sig_cent  = sig_edges(1:end-1) + diff(sig_edges)/2;
nB = numel(sig_cent);

AOUC_sig_20 = nan(nB, numel(MONTHS));
AOUC_sig_24 = nan(nB, numel(MONTHS));

for mi = 1:numel(MONTHS)
    idx20 = discretize(sig20(:,mi), sig_edges);
    tmp = nan(nB,1);
    for b = 1:nB
        k = (idx20==b);
        if any(k), tmp(b) = mean(AOU_C_molm3_20(k,mi), 'omitnan'); end
    end
    AOUC_sig_20(:,mi) = tmp;

    idx24 = discretize(sig24(:,mi), sig_edges);
    tmp = nan(nB,1);
    for b = 1:nB
        k = (idx24==b);
        if any(k), tmp(b) = mean(AOU_C_molm3_24(k,mi), 'omitnan'); end
    end
    AOUC_sig_24(:,mi) = tmp;
end

NPP20a = NPP_molC_m2_d_20 - mean(NPP_molC_m2_d_20,'omitnan');   
NPP24a = NPP_molC_m2_d_24 - mean(NPP_molC_m2_d_24,'omitnan');
AOUC_sig_20a = AOUC_sig_20 - mean(AOUC_sig_20,2,'omitnan');
AOUC_sig_24a = AOUC_sig_24 - mean(AOUC_sig_24,2,'omitnan');

lags = 0:2;
R20_sig = nan(nB, numel(lags)); P20_sig = R20_sig; N20_sig = R20_sig;
R24_sig = nan(nB, numel(lags)); P24_sig = R24_sig; N24_sig = R24_sig;

for li = 1:numel(lags)
    L = lags(li);
    x = NPP20a(1:end-L);
    for b = 1:nB
        y = AOUC_sig_20a(b,1+L:end)'; good = isfinite(x) & isfinite(y);
        N20_sig(b,li)=sum(good);
        if N20_sig(b,li)>=6, [r,p]=corr(x(good),y(good),'Type','Spearman'); R20_sig(b,li)=r; P20_sig(b,li)=p; end
    end
    x = NPP24a(1:end-L);
    for b = 1:nB
        y = AOUC_sig_24a(b,1+L:end)'; good = isfinite(x) & isfinite(y);
        N24_sig(b,li)=sum(good);
        if N24_sig(b,li)>=6, [r,p]=corr(x(good),y(good),'Type','Spearman'); R24_sig(b,li)=r; P24_sig(b,li)=p; end
    end
end

Lag_mo = lags(:);
Max_r_20 = nan(numel(lags),1); Sig_of_max_20 = Max_r_20; P_of_max_20 = Max_r_20; N_of_max_20 = Max_r_20;
Max_r_24 = nan(numel(lags),1); Sig_of_max_24 = Max_r_24; P_of_max_24 = Max_r_24; N_of_max_24 = Max_r_24;
for li=1:numel(lags)
    [mx,ib] = max(R20_sig(:,li)); if isfinite(mx)
        Max_r_20(li)=mx; Sig_of_max_20(li)=sig_cent(ib); P_of_max_20(li)=P20_sig(ib,li); N_of_max_20(li)=N20_sig(ib,li);
    end
    [mx,ib] = max(R24_sig(:,li)); if isfinite(mx)
        Max_r_24(li)=mx; Sig_of_max_24(li)=sig_cent(ib); P_of_max_24(li)=P24_sig(ib,li); N_of_max_24(li)=N24_sig(ib,li);
    end
end

T20_iso = table(Lag_mo,Sig_of_max_20,Max_r_20,P_of_max_20,N_of_max_20, ...
    'VariableNames',{'Lag_mo','Sigma0_of_max_r','Max_r','pval','Npairs'});
T24_iso = table(Lag_mo,Sig_of_max_24,Max_r_24,P_of_max_24,N_of_max_24, ...
    'VariableNames',{'Lag_mo','Sigma0_of_max_r','Max_r','pval','Npairs'});

% ---------------------- 18) Export flux from AOUc tendency -----
[~,i20L] = max(Max_r_20,[],'omitnan'); Lbest20 = Lag_mo(i20L); if ~isfinite(Lbest20), Lbest20=0; end
[~,i24L] = max(Max_r_24,[],'omitnan'); Lbest24 = Lag_mo(i24L); if ~isfinite(Lbest24), Lbest24=0; end

dA20 = AOU_C_belowEZ_m2_20(1+Lbest20:end) - AOU_C_belowEZ_m2_20(1:end-Lbest20);
export_20 = dA20 ./ max(Lbest20*dt_days,eps);             % mol C m^-2 d^-1
dA24 = AOU_C_belowEZ_m2_24(1+Lbest24:end) - AOU_C_belowEZ_m2_24(1:end-Lbest24);
export_24 = dA24 ./ max(Lbest24*dt_days,eps);

% ---------------------- 19) PEE (%) for 2020 & 2024 ------------
NPP20_align = NPP_molC_m2_d_20(1:end-Lbest20);
NPP24_align = NPP_molC_m2_d_24(1:end-Lbest24);
PEE_2020_percent = pee_percent(export_20, NPP20_align);
PEE_2024_percent = pee_percent(export_24, NPP24_align);

% Round & table
PEE_2020_percent = round(PEE_2020_percent,2);
PEE_2024_percent = round(PEE_2024_percent,2);
PEE_table = table( ["2020";"2024"], [PEE_2020_percent; PEE_2024_percent], ...
    'VariableNames', {'Year','PEE_percent'});

if save_csv
    writetable(PEE_table, fullfile(resdir,'pee_summary.csv'));
end

% ---------------------- package outputs ------------------------
arrays = struct();
arrays.zc = zc; arrays.dz = dz;
arrays.MONTHS = MONTHS;
arrays.KPAR20 = KPAR20; arrays.KPAR24 = KPAR24;
arrays.NPP_molC_m2_d_20 = NPP_molC_m2_d_20;
arrays.NPP_molC_m2_d_24 = NPP_molC_m2_d_24;
arrays.AOU_C_molm3_20   = AOU_C_molm3_20;
arrays.AOU_C_molm3_24   = AOU_C_molm3_24;
arrays.AOU_C_belowEZ_m2_20 = AOU_C_belowEZ_m2_20;
arrays.AOU_C_belowEZ_m2_24 = AOU_C_belowEZ_m2_24;

results = struct();
results.PEE_2020_percent = PEE_2020_percent;
results.PEE_2024_percent = PEE_2024_percent;
results.PEE_table        = PEE_table;
results.lag_best_2020    = Lbest20;
results.lag_best_2024    = Lbest24;
results.T20_iso          = T20_iso;
results.T24_iso          = T24_iso;
results.arrays           = arrays;

% console print
fprintf('\n=== Particle Export Efficiency (PEE) ===\n');
fprintf('2020: %.2f %% (lag=%d months)\n', PEE_2020_percent, Lbest20);
fprintf('2024: %.2f %% (lag=%d months)\n\n', PEE_2024_percent, Lbest24);

end % main function


% ========================= HELPERS ===============================
function pct = pee_percent(export_flux, npp_areal)
% 100 * mean(export)/mean(NPP), robust to NaN/zero
mexp = mean(export_flux, 'omitnan');
mnpp = mean(npp_areal,  'omitnan');
if ~isfinite(mexp) || ~isfinite(mnpp) || mnpp==0
    pct = NaN;
else
    pct = 100 * (mexp / mnpp);
end
end
