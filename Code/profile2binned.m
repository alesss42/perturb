% Bin profiles into depth bins
%
% July-2023, Pat Welch, pat@mousebrains.com

function [row, retval] = profile2binned(row, a, pars)
arguments (Input)
    row (1,:) table % row to work on
    a struct % Output of mat2profile
    pars struct % Parameters, defaults from get_info
end % arguments Input
arguments (Output)
    row table % row worked on
    retval (2,1) cell % {filename or missing, empty or binned profiles}
end % arguments Output

retval = {missing, []};

if ~row.qProfileOkay
    return;
end % if ~row.qProfileOkay

fnProf = row.fnProf;
fnBin = fullfile(pars.prof_binned_root, append(row.name, ".mat"));
row.fnBin = fnBin;

if isnewer(fnBin, fnProf)
    retval = {fnBin, []}; % We want this for combining
    fprintf("%s: %s is newer than %s\n", row.name, row.fnBin, row.fnProf);
    return;
end % if isnewer

if isempty(a)
    fprintf("Loading %s\n", row.fnProf);
    a = load(row.fnProf);
end % if isempty

method = pars.bin_method; % Which method to aggregate the data together

pInfo = a.pInfo;
profiles = a.profiles;

fprintf("%s: Binning %d profiles\n", row.name, numel(profiles));

if pars.profile_direction == "time" % Bin in time
    binSize = seconds(pars.bin_width); % Bin stepsize in (sec)
    keyName = "t";
    binFunc = @bin_by_time;
    glueFunc = @glue_lengthwise;
else % Bin by depth
    binSize = pars.bin_width; % Bin stepsize (m)
    keyName = "depth";
    binFunc = @bin_by_real;
    glueFunc = @glue_widthwise;
end % if profile_direction

casts = cell(numel(profiles),1);
for index = 1:numel(profiles)
    profile = profiles{index};

    tFast = binFunc(binSize, keyName, profile.fast, method);
    tSlow = binFunc(binSize, keyName, profile.slow, method);
    tFast = renamevars(tFast, ["n", "t", "depth"], ["n_fast", "time_fast", "depth_fast"]);
    tSlow = renamevars(tSlow, ["n", "depth"], ["n_slow", "depth_slow"]);
    casts{index} = outerjoin(tSlow, tFast, "Keys", "bin", "MergeKeys", true);
end % for index

qDrop = cellfun(@isempty, casts); % This shouldn't happend

if any(qDrop)
    casts = casts(~qDrop);
    pInfo = pInfo(~qDrop,:);
end % any qDrop

if isempty(casts)
    row.qProfileOkay = false;
    fprintf("%s: No usable casts found in %s\n", row.name, row.fnProf);
    return;
end

tbl = glueFunc("bin", casts);

binned = struct ( ...
    "tbl", tbl, ...
    "info", pInfo);

my_mk_directory(fnBin);
save(fnBin, "-struct", "binned", pars.matlab_file_format);
fprintf("%s: Saving %d profiles to %s\n", row.name, size(pInfo,1), fnBin);
retval = {fnBin, binned};
end % profile2binned