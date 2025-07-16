% code to run with everytin else

function [dInfo, tbl] = calc_chi_v2(profile, diss, pInfo, pars)
arguments (Input)
    profile struct % Profile information
    diss table % Profile informatios
    pInfo (1,:) table % Summary information about the profile
    pars struct
end % arguments Input
arguments (Output)
    dInfo (1,:) table % pInfo with extra fields
    tbl table % Tabular form of diss struct
end % arguments Output


addpath '/Users/alesanchez-rios/Documents/CRUISES/SUNRISE_PROJ/codes_to_run/odas/'
addpath '/Users/alesanchez-rios/Documents/CRUISES/SUNRISE_PROJ/codes_to_run/'
addpath '/Users/alesanchez-rios/Documents/Science_Projects/Microstructure_chi/Code_borrow/MLE/'

kappa_T = 1.4e-7; % thermal diffusivity [m^2/s]

label = sprintf("%s(%d)", pInfo.name, pInfo.index);

[dissInfo, DT_HP, AA] = mk_diss_info(profile, pars, pInfo, label);

dNames = ["speed", "T", "t", "P"]; % Names that might need trimmed

if pars.diss_trim_top && pars.trim_calculate % Trim the top of the profile
    q = dissInfo.P >= (pInfo.trim_depth + pars.diss_trim_top_offset);
    DT_HP = DT_HP(q,:);
    AA  = AA(q,:);
    for name = dNames
        dissInfo.(name) = dissInfo.(name)(q);
    end % for name
end % if trim_use

dInfo = pInfo;
tbl = table(); % Nothing calculated

diss.t = pInfo.t0 + seconds(diss.t - diss.t(1));
dissInfo.spec_length = dissInfo.diss_length;
dissInfo.gradient_method = 'high_pass';

%% Calculation chi

scalars = [DT_HP(:,1), DT_HP(:,2)]; % scalars = [profiles{Pt}.fast.T1_dT1(:)];

sp = get_scalar_spectra_odas(scalars, [], profile.fast.P_fast, datenum(profile.fast.t_fast),...
    profile.fast.speed_fast, dissInfo);

kappa = gsw_kappa(profile.slow.SA, profile.slow.theta, profile.slow.depth);

% Depth i ant th is to happen
for z_p = 1:length(sp.P)
    % zp_i =  find(sp_1.P > z_p, 1,"first");

    % finding the depth in Pats thing
    depz = sp.P(z_p);
    z_i =  find(diss.depth > depz, 1,"first");
    z_s =  find(profile.slow.depth > depz, 1,"first");

    chi_ale.chi_P(z_p) = depz;
    chi_ale.chi_lat(z_p) = dInfo.lat;
    chi_ale.chi_lon(z_p) = dInfo.lon;


    if isempty(z_i)
        z_i = length(diss.depth);
    end

    for ss = 1:2
        out = vmp_chi_MLE(sp.K(:,z_p), sp.scalar_spec(:,ss,z_p),...
            diss.K_max(z_i,ss)+10, diss.nu(z_i), kappa(z_s));

        chi_ale.chi_just(z_p, ss) = out.chi_just;

        chi_ale.chi(z_p, ss) = out.chi;
        chi_ale.e_chi(z_p, ss) = out.epsilon_chi;
        chi_ale.kB(z_p, ss) = out.kB;
        chi_ale.delta95_kB(z_p, ss) = out.delta95_kB;
        chi_ale.k_noise_cutoff(z_p, ss) = out.k_noise_cutoff;
        chi_ale.S_batcherlor(:,ss,z_p) = out.S_batchelor;
        chi_ale.k_batcherlor(:,ss,z_p) = out.k_batchelor;
    end
end



[dInfo, tbl] = mk_diss_struct_ale(pInfo, diss, chi_ale);


%% Extra functions

    function [dissInfo, DT_HP, AA] = mk_diss_info(profile, pars, pInfo, label)
        % arguments
        %     profile struct
        %     pars struct % From get_info
        %     pInfo (1,:) table
        %     label string
        % end % arguments

        fast = profile.fast; % fast variables for despiking
        fft_length_sec = pars.diss_fft_length_sec;
        diss_length_factor = pars.diss_length_fac;

        AA = table();
        for name = ["Ax", "Ay"]
            AA.(name) = my_despike(fast.(name), pInfo.fs_fast, pars, "A", ...
                sprintf("%s %s %2g", label, name, fft_length_sec), pInfo);
        end
        AA = table2array(AA);

        % Grab all the fP07 probes
        names = regexp(string(fast.Properties.VariableNames), "^gradT\d+$", "once", "match");
        names = unique(names(~ismissing(names))); % Sorted probes, assumes <10 shear probes

        DT = table(); % Space for all the shear probes
        for name = names
            DT.(name) = fast.(name);
        end % for
        DT = table2array(DT);

        DT_cut = 0.5 * 1 / fft_length_sec; % Follow Matlab manual
        [bh, ah] = butter(1, DT_cut / pInfo.fs_fast / 2, "high");
        % Do a forward filter then flip and reverse filter
        DT_HP = filter(bh, ah, DT); % Filter forwards
        DT_HP = flipud(DT_HP); % Flip forwards to backwards
        DT_HP = filter(bh, ah, DT_HP); % Filter backwards
        DT_HP = flipud(DT_HP); % Flip backwards to forward

        dissInfo = struct();
        dissInfo.fft_length = round(fft_length_sec * pInfo.fs_fast); % FFT length in bins
        dissInfo.diss_length = diss_length_factor * dissInfo.fft_length; % Dissipation length in bins

        if pars.diss_overlap_factor == 0
            dissInfo.overlap = 0; % No overlap
        else
            dissInfo.overlap = ceil(dissInfo.diss_length / pars.diss_overlap_factor); % Number of bins to overlap
        end

        for name = ["goodman", "f_limit", "fit_2_isr", "f_AA", "fit_order"]
            val = pars.(append("diss_", name));
            if isnan(val), continue; end
            dissInfo.(name) = val;
        end

        dissInfo.fs_fast = pInfo.fs_fast;
        dissInfo.fs_slow = pInfo.fs_slow;
        dissInfo.t = fast.t_fast;
        dissInfo.P = fast.P_fast;

        if ismember(pars.diss_speed_source, string(fast.Properties.VariableNames))
            dissInfo.speed = fast.(pars.diss_speed_source);
        elseif ismember(pars.diss_speed_source, string(profile.slow.Properties.VariableNames))
            dissInfo.speed = interp1(profile.slow.t, profile.slow.(pars.diss_speed_source), fast.t, "linear", "extrap");
        else
            error("diss_speed_source, %s, not in fast table", pars.diss_speed_source);
        end % if ~ismember

        if ismissing(pars.diss_T_source)
            dissInfo.T = ...
                (pars.diss_T1_norm * fast.T1_fast + pars.diss_T2_norm * fast.T2_fast) ./ ...
                (pars.diss_T1_norm + pars.diss_T2_norm);
        else % if ismissing
            TName = pars.diss_T_source; % column name for temperature source
            if ismember(TName, string(fast.Properties.VariableNames)) % A fast variable
                dissInfo.T = fast.(TName);
            elseif ismember(TName, string(profile.slow.Properties.VariableNames)) % A slow variable
                dissInfo.T = interp1(profile.slow.t, profile.slow.(TName), fast.t, "linear", "extrap");
            else
                error("diss_T_source, %s, not found in fast nor slow", TName);
            end
        end % if ismissing
    end % mk_diss_info

    function b = my_despike(a, fs, pars, codigo, tit, pInfo)
        % arguments (Input)
        %     a (:,1) {mustBeNumeric} % Vector to be despiked
        %     fs (1,1) double {mustBePositive} % Samplig frequency
        %     pars struct % Parameters, defaults from get_info
        %     codigo string % Middle field of parameter name, A or sh
        %     tit string % diagnostic title string
        %     pInfo (1,:) table % Profile summary information for this profile
        % end % arguments Input

        p = struct();
        for name = ["thresh", "smooth", "N_FS", "warning_fraction"]
            p.(name) = pars.(sprintf("despike_%s_%s", codigo, name));
        end % for name

        [b, ~, ~, raction] = despike(a, p.thresh, p.smooth, fs, round(p.N_FS * fs));

        if raction > p.warning_fraction
            fprintf("WARNING: %s spike ratio %.1f%% for profile %d in %s\n", ...
                tit, raction * 100, ...
                pInfo.index, pInfo.name);
        end % raction >
    end % my_despike

    function [dInfo, tbl] = mk_diss_struct_ale(pInfo, diss, chi_ale)
        % arguments (Input)
        %     pInfo (1,:) table % Profile information
        %     chi struct % Dissipation information from get_diss_odas
        % end % arguments Input
        % arguments(Output)
        %     dInfo (1,:) table % Dissipation input information with a table of the dissipation results, built on pInfo
        %     tbl table % Dissipation information built into a table, rows are depth/time
        % end % arguments Output

        % I don't like hardcoding names, but for a single dissipation estimate, size fails
        pNames = ["speed", "nu", "P", "T", "AOA", "epsilonMean", 'gradT1', 'gradT2'];
        npNames = ["chi_dT1_mean", "chi_dT2_mean", "chi_dT1_e1", "chi_dT1_e2", "chi_dT2_e1", "chi_dT2_e2"];
        mnpNames = "batchelor_spec";
        mnnpNames = ["tracer_spec", "tracer_spec"];
        mpNames = ["F", "K"];

        dInfo = pInfo; % We're going to add to pInfo columns with dissipation scalars
        tbl = table();
        tbl.t = diss.t; % Force as first column
        tbl.depth = diss.depth; % Force to the second column

        for name = string(fieldnames(diss))'
            
            
            if ismember(name, pNames)
                 val = diss.(name);
                
                tbl.(name) = val;
            elseif name == 'Properties'
            elseif name == 'Row'
            elseif name == 'Variables'
            else
                
            end
        end



            tbl.chi_dT1_mean =  nanmean(chi_ale.chi_just,2);
            tbl.chi_dT2_mean =  nanmean(chi_ale.chi,2);
            tbl.chi_dT1_e1 = chi_ale.chi_just(:,1);
            tbl.chi_dT1_e2 = chi_ale.chi_just(:,2);
            tbl.chi_dT2_e1 = chi_ale.chi(:,1);
            tbl.chi_dT2_e2 = chi_ale.chi(:,2);
            tbl.chi_e1 = chi_ale.e_chi(:,1);
            tbl.chi_e2 = chi_ale.e_chi(:,2);
            tbl.bathchelor_spec = permute(chi_ale.S_batcherlor, [3, 2, 1]);
            tbl.bathchelor_k = permute(chi_ale.k_batcherlor, [3, 2, 1]);
            tbl.scalar_spec = permute(sp.scalar_spec, [3,2,1]);
            tbl.k_noice_cutoff = chi_ale.k_noise_cutoff;
            tlb.kB = chi_ale.kB;

        
    end % mk_diss_struct
end