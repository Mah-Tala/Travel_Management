clc
clear all
close all

% -------------------------------------------------------------------------
% generation_tau_MA.m
% Purpose: build time_outside and tau for MA counties from daily flows.
% Inputs (relative to this script's folder Code/data_generation):
%   ../../Datasets/cbg_fips_codes.csv
%   ../../Datasets/Massachusetts_county/April/population.mat   (variable: population)
%   ../../Datasets/2020_flows/daily_county2county_2020_08_01.csv
% Output (DO NOT CHANGE NAME):
%   ../../Datasets/Massachusetts_county/April/travel.mat
%     variables: time_outside, MA_county_id, tau
% -------------------------------------------------------------------------

% --- County IDs for Massachusetts (state_fips == 25) ---------------------
cbg = '../../Datasets/cbg_fips_codes.csv';
T_cbg = readtable(cbg);
county_id = double(T_cbg.county_fips);
state_id  = double(T_cbg.state_fips);

MA_county_id = county_id(find(state_id == 25));   % MA county FIPS (as used downstream)
county_num   = length(MA_county_id);

% --- Population vector (aligned with your MA_county_id order) -------------
load('../../Datasets/Massachusetts_county/April/population.mat');  % loads 'population'

% --- Daily county-to-county flows (pick the date file you want) -----------
flows = readtable('../../Datasets/2020_flows/county2county/daily_county2county_2020_08_01.csv');

pop_flows     = zeros(county_num, county_num);
visitor_flows = zeros(county_num, county_num);

origin_id_all      = double(flows.geoid_o);
destination_id_all = double(flows.geoid_d);
flows_all          = double(flows.pop_flows);
visits             = double(flows.visitor_flows);

% Keep rows where the origin is in MA (state prefix 25xxx)
MA_origin_index = find(floor(origin_id_all/1e3) == 25);
MA_origins      = origin_id_all(MA_origin_index);
MA_destinations = destination_id_all(MA_origin_index);
MA_flows        = flows_all(MA_origin_index);
MA_visits       = visits(MA_origin_index);

% --- Build origin->destination matrices for MA counties -------------------
for i = 1:county_num
    idx_o = find(MA_origins == (25e3 + MA_county_id(i)));
    if (size(idx_o,1) ~= 0)
        flows_p = MA_flows(idx_o);
        flows_v = MA_visits(idx_o);
        for j = 1:county_num
            idx_d = find(MA_destinations(idx_o) == (25e3 + MA_county_id(j)));
            if (size(idx_d,1) ~= 0)
                pop_flows(i,j)     = flows_p(idx_d);
                visitor_flows(i,j) = flows_v(idx_d);
            end
        end
    end
end

% --- Convert flows to time spent outside (assumptions below) --------------
in_county_trip   = 1;  % time weight for trips within the same county
out_conty_traip  = 3;  % (typo kept as-is) time weight for trips to other counties

diag_inds    = logical(eye(size(pop_flows)));  % diagonal (i==j)
offdiag_inds = ~diag_inds;                     % off-diagonal (i~=j)

time_outside = pop_flows;
time_outside(diag_inds)    = pop_flows(diag_inds)    * in_county_trip;
time_outside(offdiag_inds) = pop_flows(offdiag_inds) * out_conty_traip;

% --- Tau: time outside per capita ----------------------------------------
tau = time_outside ./ population;

% --- Save (names unchanged for pipeline compatibility) -------------------
save_name = '../../Datasets/Massachusetts_county/April/travel.mat';
save(save_name, 'time_outside', 'MA_county_id', 'tau');
