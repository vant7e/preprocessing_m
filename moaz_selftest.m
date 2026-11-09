function moaz_selftest(letswave_root)
% Run on the compute node before the real BDF. Does not fit ICA.
addpath(genpath(letswave_root));
h=struct('filetype','time_amplitude','name','selftest','tags','','history',[], ...
    'datasize',[1 2 1 1 1 2048],'xstart',0,'xstep',1/512, ...
    'ystart',0,'ystep',1,'zstart',0,'zstep',1);
h.chanlocs=struct('labels',{'Fp1','Fp2'},'topo_enabled',0,'SEEG_enabled',0);
h.events=struct('code',{'1','2'},'latency',{1,2},'epoch',{1,1});
t=(0:2047)/512; d=zeros(h.datasize); d(1,1,1,1,1,:)=sin(2*pi*10*t);
d(1,2,1,1,1,:)=cos(2*pi*12*t);
[eh,ed]=RLW_segmentation(h,d,{'1','2'},'x_start',-0.1,'x_duration',1);
assert(isequal(eh.datasize,[2 2 1 1 1 513]));
first=round((1-0.1)*512)+1;
assert(isequal(squeeze(ed(1,1,1,1,1,:)),squeeze(d(1,1,1,1,1,first:first+512))));
[~,dc]=RLW_dc_removal(eh,ed+12,'linear_detrend',0);
a=mean(dc,6); assert(max(abs(a(:)))<1e-10,'DC removal mismatch.');
% +/-700 is retained; 700.01 is rejected. This is not a peak-to-peak limit.
x=zeros(eh.datasize); x(1,1,1,1,1,1:2)=[700 -700];
x(2,1,1,1,1,1:2)=[700.01 -700.01];
[ah,ad]=RLW_reject_epochs_amplitude(eh,x,'criterion',700);
assert(size(ad,1)==1&&isequal(ad,x(1,:,:,:,:,:))); assert(ah.datasize(1)==1);
[~,fd]=RLW_butterworth_filter(h,d,'filter_type','bandpass','low_cutoff',0.1,'high_cutoff',40,'filter_order',4);
[b,a]=butter(2,[0.1 40]/256,'bandpass'); ref=filtfilt(b,a,sin(2*pi*10*t));
err=squeeze(fd(1,1,1,1,1,:))-ref(:); assert(max(abs(err))<1e-10,'Native filter order semantics mismatch.');
% One valid + one out-of-bounds event: native's zero tail must be removed.
h.events(2).latency=3.5;
[bh,bd,valid]=moaz_segment_trials(h,d,{'1','2'},'drop_incomplete_with_audit');
assert(isequal(valid,[true false])&&size(bd,1)==1&&bh.datasize(1)==1);
assert(isequal(bd,ed(1,:,:,:,:,:)));
assert(round(291.1*512)+1==149044);
c=LW_ICA_compute_merged('default');
assert(strcmp(c.parameters.ICA_algorithm,'runica')&&strcmp(c.parameters.ICA_mode,'square'));
root=tempname; mkdir(root); cleaner=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
CLW_save(root,eh,ed); saved=load(fullfile(root,'selftest.mat'),'data');
assert(isa(saved.data,'single')&&isequal(saved.data,single(ed)));
fprintf('Native Moaz preflight self-test PASSED. No real EEG processed.\n');
end
