clc
clear all

% -------------------------------------------------------------------------
% generation_initial_rate_MA.m
% Purpose:
%   Build April 1, 2020 initial rates for MA counties:
%   s_init (susceptible), h_init (healed), x_init (infected), d_init (deceased)
% Inputs (relative to this script's folder Code/data_generation):
%   ../../Datasets/us-counties-2020.csv
%   ../../Datasets/cbg_fips_codes.csv
%   ../../Datasets/us-counties_1.csv
%   ../../Datasets/Massachusetts_county/April/population.mat   (vector "population")
% Output (DO NOT CHANGE NAME):
%   ../../Datasets/Massachusetts_county/April/initial_rate_04_01.mat
%     variables: s_init, h_init, x_init, d_init, cum_num_county, recovery_num_county, death_num_county
% -------------------------------------------------------------------------

%% ------------------------- Cases on 2020-04-01 ---------------------------
% Read cumulative COVID cases and filter to Massachusetts on 01-Apr-2020
name = '../../Datasets/us-counties-2020.csv';
T_test = readtable(name);

test_date_all     = string(T_test.date);
cum_num_all       = double(T_test.cases);
county_names_all_1 = string(T_test.county);
states_all        = string(T_test.state);

% Get MA rows on 01-Apr-2020
index = find(test_date_all == '01-Apr-2020' & states_all == 'Massachusetts');
test_date    = test_date_all(index);   % not used later, kept for clarity
cum_num      = cum_num_all(index);
county_names_1 = county_names_all_1(index);

% Match naming style to FIPS (append " County")
for i = 1:length(county_names_1)
    county_names_1(i) = join([county_names_1(i), "County"], " ");
end

%% ---------------------- FIPS county order for MA -------------------------
% Use FIPS to get the canonical county order and names (state_fips == 25)
name2 = '../../Datasets/cbg_fips_codes.csv';
T_fips_code = readtable(name2);

state_ids         = T_fips_code.state_fips;
county_ids        = T_fips_code.county_fips;
county_names_all_2 = string(T_fips_code.county);

index          = find(state_ids == 25);
MA_county_ids  = county_ids(index);          % keep if used downstream
county_names_2 = county_names_all_2(index);

n = length(county_names_2);  % number of MA counties

% Map cumulative cases to the FIPS county order
cum_num_county = zeros(n, 1);
for i = 1:n
    county_index = find(county_names_1 == county_names_2(i));
    cum_num_county(i) = cum_num(county_index);  % cumulative confirmed per county
end

%% -------------------------- Load population ------------------------------
% 'population' must be aligned to the same MA county order as county_names_2
name3 = '../../Datasets/Massachusetts_county/April/population.mat';
load(name3)  % loads variable: population

% Estimated reporting rate = 0.14  → scale reported counts to total
s_0 = 1 - (cum_num_county / 0.14) ./ population;  % susceptible rate per county

clear test_date_all cum_num_all county_name_all_1 index  % lightweight cleanup

%% ---------------------------- Deaths (Apr 1) -----------------------------
% Read deaths and filter to MA on 01-Apr-2020
name = '../../Datasets/us-counties-2020.csv';
T_death = readtable(name);

death_date_all  = string(T_death.date);
death_num_all   = double(T_death.deaths);
county_names_all = string(T_death.county);
state_names_all = string(T_death.state);

index_date      = find(death_date_all == '01-Apr-2020');   % Apr 1, 2020
death_num_date  = death_num_all(index_date);
state_name_date = state_names_all(index_date);
county_names_date = county_names_all(index_date);

index_state = find(state_name_date == 'Massachusetts');    % MA only
death_num   = death_num_date(index_state);
county_names = county_names_date(index_state);

% Match naming style to FIPS (append " County")
for i = 1:length(county_names)
    county_names(i) = join([county_names(i), "County"], " ");
end

% Map deaths to the FIPS county order; if a county missing → keep Inf
death_num_county = inf * ones(n, 1);
for i = 1:n
    index_death = find(county_names == county_names_2(i));
    if ~isempty(index_death)
        death_num_county(i) = death_num(index_death);
    end
end

clear county_names

%% ----------------------- Recovery and initial rates ----------------------
% Rough recovery estimate based on ratio 8878/215215 (kept as in original)
recovery_num_county = round(cum_num_county * 8878 / 215215);

% Reporting rate again (0.14) to scale to totals; split infectious vs deaths
s_init = s_0;
h_init = ((death_num_county + recovery_num_county) / 0.14) ./ population;                 % healed
x_init = ( (cum_num_county - death_num_county - recovery_num_county) / 0.14 * 0.86 ) ./ population; % infectious
d_init = ( (cum_num_county - death_num_county - recovery_num_county) / 0.14 * 0.14 ) ./ population; % deceased

% Quick sanity print (sum of compartments per county)
s_init + h_init + x_init + d_init;

% ----------------------------- Save outputs ------------------------------
save_name = ['../../Datasets/Massachusetts_county/April/initial_rate_04_01', '.mat'];
save(save_name, 's_init', 'h_init', 'x_init', 'd_init', 'cum_num_county', 'recovery_num_county', 'death_num_county');
