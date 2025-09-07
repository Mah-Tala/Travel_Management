% Build MA county population vector (aligned to FIPS county order) -> population.mat
clc; clear; close all;

% 1) Read population Excel (columns: "Geographic Area", "2020 Census")
excelFile = '../../Datasets/Population_MA.xlsx';
popTable  = readtable(excelFile, 'VariableNamingRule','preserve');

countyNamesPop = string(popTable.("Geographic Area"));
pop2020        = popTable.("2020 Census");
if ~isnumeric(pop2020)
    pop2020 = str2double(replace(string(pop2020), ",", "")); % remove commas then to double
end

% 2) Read FIPS and get MA (state_fips == 25) county names in the desired order
fipsFile  = '../../Datasets/cbg_fips_codes.csv';
fipsTable = readtable(fipsFile, 'VariableNamingRule','preserve');

isMA            = (fipsTable.state_fips == 25);
countyIdsMA     = fipsTable.county_fips(isMA);
countyNamesFips = string(fipsTable.county(isMA));

% keep unique counties while preserving FIPS order (cbg file may repeat counties)
[countyIdsMA, ia] = unique(countyIdsMA, 'stable');
countyNamesFips    = countyNamesFips(ia);

% 3) Map population to the FIPS county order and save the vector
% (light name cleaning to ensure matches)
cnPopClean  = lower(strtrim(countyNamesPop));
cnPopClean  = replace(cnPopClean, " county", "");
cnFipsClean = lower(strtrim(countyNamesFips));
cnFipsClean = replace(cnFipsClean, " county", "");

population = zeros(numel(cnFipsClean), 1);
for i = 1:numel(cnFipsClean)
    j = find(cnPopClean == cnFipsClean(i), 1);
    population(i) = pop2020(j);
end

outDir = '../../Datasets/Massachusetts_county/April';
if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, 'population.mat'), 'population');
