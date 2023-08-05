%
% Convert test .P files to binned data
%
% August-2023, Pat Welch, pat@mousebrains.com

my_root = fileparts(mfilename("fullpath"));
code_root = fullfile(my_root, "../Code");
data_root = fullfile(my_root, "../Data");
vmp_root = fullfile(data_root, "VMP");
output_root = fullfile(my_root, "../Temp");
GPS_filename = fullfile(data_root, "GPS/gps.csv");

addpath(code_root, "-begin");

process_VMP_files( ...
    "debug", true, ...
    "vmp_root", vmp_root, ...      % Where the input .P files are located
    "output_root", output_root, ... % Where to write output to
    "gps_class", @GPS_from_csv, ... % Class to supply GPS data
    "gps_filename", GPS_filename, ...  % Input GPS data
    "diss_epsilon_minimum", 1e-11, ...   % Drop dissipation estimates smaller than this value
    "netCDF_contributor_name", "Pat Welch", ...
    "netCDF_contributor_role", "researcher", ...
    "netCDF_creator_name", "Pat Welch", ...
    "netCDF_creator_email", "pat@mousebrains.com", ...
    "netCDF_creator_institution", "CEOAS, Oregon State University", ...
    "netCDF_creator_type", "researcher", ...
    "netCDF_creator_url", "https://arcterx.ceoas.oregonstat.edu", ...
    "netCDF_id", "Test of software suite", ...
    "netCDF_institution", "CEOAS, Oregon State University", ...
    "netCDF_platform", "Rockland VMP250", ...
    "netCDF_product_version", "0.1", ...
    "netCDF_program", "Test of software suite", ...
    "netCDF_project", "Test of software suite" ...
    );

rmpath(code_root);