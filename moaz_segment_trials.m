function [out_header,out_data,valid]=moaz_segment_trials(header,data,event_labels,policy)
% Explicit compatibility repair: trim native preallocated zero tail epochs.
% Preserve the fixed block duration and all native event/rounding conventions.
idx=[];
for k=1:numel(header.events)
    if any(strcmpi(header.events(k).code,event_labels)), idx(end+1)=k; end %#ok<AGROW>
end
valid=false(size(idx)); n=round(1/header.xstep)+1;
for k=1:numel(idx)
    first=round((header.events(idx(k)).latency-0.1-header.xstart)/header.xstep)+1;
    valid(k)=first>=1&&first+n-1<=header.datasize(6);
end
if any(~valid)
    assert(strcmp(policy,'drop_incomplete_with_audit'), ...
        'Incomplete trials at block boundary. Explicit boundary policy required.');
end
assert(any(valid),'No complete trials remain in this block.');
[out_header,out_data]=RLW_segmentation(header,data,event_labels,'x_start',-0.1,'x_duration',1);
% This pinned native implementation packs valid trials first, but leaves a
% preallocated zero tail for invalid trials. Never let that tail reach ICA.
assert(size(out_data,1)==numel(idx),'Unexpected native segmentation allocation.');
if any(~valid)
    tail=out_data(sum(valid)+1:end,:,:,:,:,:);
    assert(all(tail(:)==0),'Native zero-tail behavior changed; refusing repair.');
end
out_data=out_data(1:sum(valid),:,:,:,:,:);
out_header.datasize(1)=sum(valid);
if ~isempty(out_header.events)
    assert(all([out_header.events.epoch]>=1&[out_header.events.epoch]<=sum(valid)), ...
        'Unexpected native event epochs after boundary repair.');
end
end
