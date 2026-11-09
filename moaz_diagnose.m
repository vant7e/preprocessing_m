function moaz_diagnose(config_path)
% Native Letswave 6 protocol branch. Deliberately stops before ICA application.
cfg=jsondecode(fileread(config_path));
if ~isfield(cfg,'boundary_policy'), cfg.boundary_policy='error'; end
assert(~isfolder(cfg.output)&&~isfile(cfg.output),'Output already exists; choose a NEW output directory.');
assert(isfile(cfg.bdf),'BDF does not exist: %s',cfg.bdf);
assert(strcmp(cfg.subject,'subj04')&&strcmp(cfg.session,'S1'),'This pilot configuration is restricted to subj04 S1.');
mkdir(cfg.output); mkdir(fullfile(cfg.output,'qc')); mkdir(fullfile(cfg.output,'datasets')); mkdir(fullfile(cfg.output,'ica_diagnostics'));
diary(fullfile(cfg.output,'matlab.log')); cleanup=onCleanup(@() diary('off')); %#ok<NASGU>
manifest=struct('status','running','subject',cfg.subject,'session',cfg.session,'configuration',cfg, ...
    'matlab_version',version,'started',datestr(now,30),'ica_applied',false,'excluded',[], ...
    'historical_bitwise_equivalence',false);
write_json(fullfile(cfg.output,'run.json'),manifest);
try
    addpath(genpath(cfg.letswave_root));
    here=fileparts(mfilename('fullpath'));
    lock=jsondecode(fileread(fullfile(here,'source_lock.json')));
    for k=1:numel(lock.files)
        f=lock.files(k); actual=which(f.function);
        expected=fullfile(cfg.letswave_root,f.relative_path);
        assert(~isempty(actual),'Missing native function: %s',f.function);
        assert(strcmp(canonical(actual),canonical(expected)),'Wrong/shadowed Letswave function %s: %s',f.function,actual);
        assert(strcmp(hash_file(actual),f.sha256),'Source mismatch: %s. Use the pinned checkout.',actual);
    end
    dependencies=struct();
    for f={'readlocs','topoplot','ft_read_data','ft_read_header','ft_read_event','butter','filtfilt'}
        actual=which(f{1});
        assert(~isempty(actual),'Missing dependency: %s',f{1});
        dependencies.(f{1})=struct('path',actual,'sha256',hash_file(actual));
        if ~ismember(f{1},{'butter','filtfilt'})
            prefix=[canonical(cfg.letswave_root) filesep];
            assert(startsWith(canonical(actual),prefix),'Shadowed native dependency: %s',actual);
        end
    end
    moaz_selftest(cfg.letswave_root);
    write_json(fullfile(cfg.output,'source_lock_used.json'),lock);
    copyfile(config_path,fullfile(cfg.output,'config_used.json'));
    manifest.source_commit=lock.commit;
    manifest.dependencies=dependencies;
    manifest.wrapper_sha256=hash_file([mfilename('fullpath') '.m']);
    manifest.limitations={ ...
        'Moaz historical Letswave version and Butterworth order are unknown; pinned LW6 default order 4 is an assumption.', ...
        'Native CLW_save single-precision round trips are retained between stages.', ...
        'Native RLW_ICA_compute prepends a zero array equal in length to the real data; unmodified for this source version.', ...
        'Native runica resets its random state from the clock; independent refits are not bitwise reproducible.', ...
        '291.1 seconds produces 149044 samples at 512 Hz in this native version; Moaz also wrote 149043.', ...
        'Native BDF event sample-to-latency convention is retained without a one-sample correction.', ...
        'No manual visual rejection/interpolation/extra ghost-event deletion has been inferred.', ...
        'Incomplete boundary trials are excluded only under the explicit drop_incomplete_with_audit policy; native zero tail is trimmed and exclusions logged.', ...
        'Diagnostics only: no ICA application, final detrend, rereference, baseline or ERP export yet.'};
    write_json(fullfile(cfg.output,'run.json'),manifest);
    
    fprintf('Importing original BDF with native Letswave...\n');
    [h,d]=RLW_import_BDF(cfg.bdf);
    assert(h.datasize(1)==1,'Expected one continuous recording.');
    assert(abs(1/h.xstep-512)<1e-9,'Expected 512 Hz; no resampling is allowed in this branch.');
    assert(h.xstart==0,'Unexpected continuous xstart; audit timing before continuing.');
    [h,d]=stage(cfg,'01_import',h,d);
    labels={h.chanlocs.labels};
    remove=strcmpi(labels,'Status');
    for k=1:numel(labels)
        remove(k)=remove(k)||~isempty(regexpi(strtrim(labels{k}),'^(EX|EXG)[1-8]$','once'));
    end
    keep=find(~remove);
    assert(numel(keep)==64,'Expected 64 scalp channels after EX/EXG/Status removal; found %d. Check labels.',numel(keep));
    assert(numel(unique(lower(string(labels(keep)))))==64,'Duplicate channel labels.');
    verify_bdf_units(cfg.bdf,labels(keep));
    h.chanlocs=h.chanlocs(keep); d=d(:,keep,:,:,:,:); h.datasize(2)=numel(keep);
    locfile=fullfile(cfg.letswave_root,'core_functions','Standard-10-20-Cap81.locs');
    locs=readlocs(locfile);
    for k=1:numel(locs), locs(k).topo_enabled=1; locs(k).SEEG_enabled=0; end
    [h,~]=RLW_edit_electrodes(h,locs);
    for k=1:64
        assert(isfield(h.chanlocs(k),'X')&&~isempty(h.chanlocs(k).X), ...
            'No native 10/20 template coordinate for %s',h.chanlocs(k).labels);
    end
    [h,d]=stage(cfg,'02_channels_coordinates',h,d);
    
    original=h.events; codes=event_codes(original);
    original_events=table((1:numel(original))',string({original.code})', ...
        [original.latency]',round([original.latency]'/h.xstep), ...
        'VariableNames',{'native_event_index','code','native_latency_s','native_latency_div_xstep'});
    writetable(original_events,fullfile(cfg.output,'qc','events_original.csv'));
    del=ismember(codes,cfg.remove_event_codes);
    % No heuristic event repair: any future manual deletion must be explicitly configured.
    assert(isempty(cfg.manual_delete_event_indices),'Pilot manual deletion list must remain empty until reviewed.');
    deleted=original_events(del,:); writetable(deleted,fullfile(cfg.output,'qc','events_deleted.csv'));
    h.events=original(~del);
    [h,d]=stage(cfg,'03_events_cleaned',h,d);
    codes=event_codes(h.events);
    triggers=arrayfun(@num2str,202:217,'UniformOutput',false);
    for k=1:16, assert(sum(codes==201+k)==1,'Expected exactly one block trigger %d.',201+k); end
    selected=ismember(codes,[198 199 201:217]);
    trigger_table=table(codes(selected)',[h.events(selected).latency]', ...
        'VariableNames',{'code','native_latency_s'});
    writetable(trigger_table,fullfile(cfg.output,'qc','block_triggers.csv'));
    fprintf('Filtering ENTIRE continuous recording 0.1-40 Hz, native order 4...\n');
    [h,d]=RLW_butterworth_filter(h,d,'filter_type','bandpass', ...
        'low_cutoff',0.1,'high_cutoff',40,'filter_order',4);
    assert(all(isfinite(d(:))),'Nonfinite filtered samples.');
    [h,d]=stage(cfg,'04_butterworth',h,d);
    
    datasets=struct('header',{},'data',{}); ledger=table(); block_summary=table(); identity_counts=table(); boundary_excluded=table();
    faces=arrayfun(@num2str,1:60,'UniformOutput',false);
    for b=1:16
        code=201+b; anchor=find(codes==code);
        block_latency=h.events(anchor).latency;
        check_bounds(h,anchor,-0.02,cfg.block_duration_s);
        [bh,bd]=RLW_segmentation(h,d,triggers(b),'x_start',-0.02,'x_duration',cfg.block_duration_s);
        [bh,bd]=stage(cfg,sprintf('05_blk_%d',code),bh,bd);
        bc=event_codes(bh.events); fi=find(ismember(bc,1:60));
        assert(~isempty(fi),'Block %d contains no ordinary faces.',code);
        all_fi=fi;
        [th,td,valid]=moaz_segment_trials(bh,bd,faces,cfg.boundary_policy);
        invalid_fi=all_fi(~valid);
        excluded_rows=table(repmat(code,numel(invalid_fi),1),bc(invalid_fi)', ...
            block_latency+reshape([bh.events(invalid_fi).latency],[],1), ...
            repmat("incomplete_trial_at_fixed_block_boundary",numel(invalid_fi),1), ...
            'VariableNames',{'block_code','identity','native_stimulus_latency_s','reason'});
        boundary_excluded=[boundary_excluded;excluded_rows]; %#ok<AGROW>
        writetable(boundary_excluded,fullfile(cfg.output,'qc','boundary_excluded_trials.csv'));
        fi=all_fi(valid);
        assert(th.datasize(1)==numel(fi)&&th.datasize(6)==513,'Unexpected native epoch shape.');
        [th,td]=stage(cfg,sprintf('06_trl_%d',code),th,td);
        [th,td]=RLW_dc_removal(th,td,'linear_detrend',0);
        [th,td]=stage(cfg,sprintf('07_dc_%d',code),th,td);
        peak=max(abs(reshape(td,size(td,1),[])),[],2);
        accepted=peak<=700;
        race=repmat("Asian",numel(fi),1); race(bc(fi)>30)="Caucasian";
        rows=table(repmat(code,numel(fi),1),(1:numel(fi))',bc(fi)',race, ...
            block_latency+[bh.events(fi).latency]',double(peak),accepted, ...
            'VariableNames',{'block_code','trial_in_block','identity','stimulus_race', ...
            'native_stimulus_latency_s','peak_abs_after_dc_uv','accepted_pre_ica'});
        ledger=[ledger;rows]; %#ok<AGROW>
        [ah,ad]=RLW_reject_epochs_amplitude(th,td,'criterion',700,'select_channels',0, ...
            'x_limits',0,'y_limits',0,'z_limits',0);
        assert(size(ad,1)==sum(accepted),'Native amplitude rejection differs from audit.');
        assert(size(ad,1)>0,'No retained trials in block %d.',code);
        [ah,ad]=stage(cfg,sprintf('08_ar_%d',code),ah,ad);
        datasets(b).header=ah; datasets(b).data=ad;
        block_summary=[block_summary; table(code,numel(all_fi),sum(accepted),sum(~accepted),sum(~valid), ...
            sum(ismember(bc,101:160)),numel(bd), ...
            'VariableNames',{'block_code','ordinary_before','ordinary_kept','ordinary_rejected_amplitude','ordinary_rejected_boundary', ...
            'oddballs_inside_block','block_data_elements'})]; %#ok<AGROW>
        for id=1:60
            identity_counts=[identity_counts;table(code,id,sum(bc(all_fi)==id),sum(accepted&(bc(fi)'==id)),sum(bc(invalid_fi)==id), ...
                'VariableNames',{'block_code','identity','n_before','n_kept','n_boundary_excluded'})]; %#ok<AGROW>
        end
        fprintf('Block %d: %d -> %d trials (DC then |amplitude| > 700 uV).\n',code,numel(fi),sum(accepted));
        fprintf('  Boundary exclusions: %d; amplitude exclusions: %d.\n',sum(~valid),sum(~accepted));
        % Write incrementally so an error never conceals the earlier audit.
        writetable(ledger,fullfile(cfg.output,'qc','trial_ledger.csv'));
        writetable(block_summary,fullfile(cfg.output,'qc','block_counts.csv'));
        writetable(identity_counts,fullfile(cfg.output,'qc','identity_block_counts.csv'));
    end
    clear d bd td ad
    if any(identity_counts.n_before~=4)
        warning('Some identity/block counts differ from four. Inspect identity_block_counts.csv; no events have been fabricated or dropped to force four.');
    end
    cfg_ica=LW_ICA_compute_merged('default');
    save(fullfile(cfg.output,'pre_ica_cache.mat'),'datasets','cfg_ica','cfg','-v7.3');
    fprintf('Fitting native merged-session ICA; no components will be excluded.\n');
    rng_before=rng; started=tic;
    [~,ica_datasets]=LW_ICA_compute_merged('process',cfg_ica,datasets,[]);
    elapsed=toc(started); rng_after=rng;
    ica_configuration=ica_datasets(1).header.history(end).configuration;
    ICA_mm=ica_configuration.parameters.ICA_mm; ICA_um=ica_configuration.parameters.ICA_um;
    assert(all(isfinite(ICA_mm(:)))&&all(isfinite(ICA_um(:))),'Nonfinite ICA matrices.');
    save(fullfile(cfg.output,'ica_matrix.mat'),'ICA_mm','ICA_um','ica_configuration', ...
        'rng_before','rng_after','elapsed','-v7.3');
    for b=1:16
        [~,~]=stage(cfg,sprintf('09_ica_%d',201+b),ica_datasets(b).header,ica_datasets(b).data);
    end
    manifest.status='ica_fitted_diagnostics_pending';
    manifest.n_trials_boundary_excluded=height(boundary_excluded);
    manifest.n_trials_before=height(ledger)+height(boundary_excluded); manifest.n_trials_kept=sum(ledger.accepted_pre_ica);
    manifest.ica_fit_seconds=elapsed; manifest.n_components=size(ICA_um,1);
    write_json(fullfile(cfg.output,'run.json'),manifest);
    moaz_plot_ica(cfg.output,datasets,ICA_mm,ICA_um);
    review=struct('subject',cfg.subject,'session',cfg.session,'source_output',cfg.output, ...
        'component_indexing','Letswave ONE-BASED: ICA1 is component 1, not Python index 1', ...
        'review_confirmed',false,'excluded_components_one_based',[], ...
        'eog_confirmed',false,'reviewer','','notes','No automatic rejection. Frontal correlation is diagnostic only.');
    write_json(fullfile(cfg.output,'ica_review_template.json'),review);
    manifest.status='awaiting_ica_review'; manifest.finished=datestr(now,30);
    write_json(fullfile(cfg.output,'run.json'),manifest);
    fprintf('DONE: awaiting manual ICA review. No final cleaned epochs have been exported.\n%s\n',cfg.output);
catch ME
    manifest.status='failed'; manifest.error=getReport(ME,'extended','hyperlinks','off');
    write_json(fullfile(cfg.output,'run.json'),manifest); rethrow(ME);
end
end

function [h,d]=stage(cfg,name,h,d)
% Match native GUI storage precision and native header sanitation between steps.
h.name=name; root=fullfile(cfg.output,'datasets');
assert(~isfile(fullfile(root,[name '.lw6']))&&~isfile(fullfile(root,[name '.mat'])),'Stage already exists.');
CLW_save(root,h,d);
a=load(fullfile(root,[name '.lw6']),'-mat'); h=a.header;
a=load(fullfile(root,[name '.mat']),'data'); d=a.data;
assert(all(isfinite(d(:))),'Nonfinite samples in %s',name);
end

function c=event_codes(events)
c=nan(1,numel(events));
for i=1:numel(events), c(i)=str2double(events(i).code); end
end

function check_bounds(h,idx,start,duration)
n=round(duration/h.xstep)+1;
for k=idx(:)'
    first=round((h.events(k).latency+start-h.xstart)/h.xstep)+1;
    assert(first>=1&&first+n-1<=h.datasize(6), ...
        'Native segmentation would exceed data: event %s latency %.9f. Audit before proceeding.', ...
        h.events(k).code,h.events(k).latency);
end
end

function verify_bdf_units(path,labels)
% Confirm threshold units from the original EDF/BDF header; do not guess.
fid=fopen(path,'r'); assert(fid>=0); close_file=onCleanup(@() fclose(fid)); %#ok<NASGU>
fseek(fid,252,'bof'); n=str2double(char(fread(fid,4,'*char')'));
assert(isfinite(n)&&n>0,'Invalid BDF channel header.');
fseek(fid,256,'bof'); names=cellstr(reshape(char(fread(fid,16*n,'*char')),16,n)');
fseek(fid,256+n*(16+80),'bof'); units=cellstr(reshape(char(fread(fid,8*n,'*char')),8,n)');
for k=1:numel(labels)
    idx=find(strcmp(strtrim(names),strtrim(labels{k}))); assert(numel(idx)==1,'BDF label mismatch.');
    unit=lower(strtrim(units{idx}));
    assert(any(strcmp(unit,{'uv','microv'})),'Expected BDF physical units uV for %s; got %s. No silent rescaling.',labels{k},unit);
end
end

function value=hash_file(path)
fid=fopen(path,'r'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8'); digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes); value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end

function write_json(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end

function p=canonical(path)
p=char(java.io.File(path).getCanonicalPath());
end
