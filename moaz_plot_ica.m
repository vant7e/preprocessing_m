function moaz_plot_ica(output,datasets,ICA_mm,ICA_um)
% Diagnostic-only. Index labels are native Letswave ONE-BASED.
root=fullfile(output,'ica_diagnostics'); if ~isfolder(root), mkdir(root); end
h=datasets(1).header; labels={h.chanlocs.labels}; nc=size(ICA_um,1);
d=cat(1,datasets.data); nt=size(d,6); ne=size(d,1);
X=double(reshape(permute(d,[2 6 1 3 4 5]),size(d,2),[])); clear d
S=ICA_um*X;
fp1=find(strcmpi(labels,'Fp1')); fp2=find(strcmpi(labels,'Fp2'));
assert(numel(fp1)==1&&numel(fp2)==1,'Fp1/Fp2 diagnostic reference missing.');
ref1=X(fp1,:); ref2=X(fp2,:); clear X
ref1=ref1-mean(ref1); ref2=ref2-mean(ref2);
metrics=zeros(nc,4); time=h.xstart+(0:nt-1)*h.xstep;
freq=(0:floor(nt/2))/(nt*h.xstep); window=0.5-0.5*cos(2*pi*(0:nt-1)/(nt-1));
for ic=1:nc
    fprintf('ICA diagnostics %d/%d\n',ic,nc);
    s=S(ic,:); centered=s-mean(s);
    r1=(centered*ref1')/max(eps,norm(centered)*norm(ref1));
    r2=(centered*ref2')/max(eps,norm(centered)*norm(ref2));
    metrics(ic,:)=[var(s),r1,r2,max(abs([r1 r2]))];
    trials=reshape(s,nt,ne)';
    fig=figure('Visible','off','Color','w','Position',[50 50 1100 780]);
    subplot(2,2,1);
    topoplot(ICA_mm(:,ic),h.chanlocs,'electrodes','on');
    title(sprintf('ICA%d | native one-based',ic)); colorbar;
    subplot(2,2,2);
    imagesc(time,1:ne,trials); axis xy; xlabel('Time from face onset (s)'); ylabel('Retained trial');
    vals=sort(abs(s)); lim=vals(max(1,ceil(0.99*numel(vals))));
    if lim>0, caxis([-lim lim]); end
    title('All retained trials; display clipped at 99th percentile'); colorbar;
    subplot(2,2,3);
    centered_trials=trials-mean(trials,2);
    spectrum=fft(centered_trials.*window,[],2);
    power=mean(abs(spectrum(:,1:numel(freq))).^2,1)/(sum(window.^2)/h.xstep);
    plot(freq,10*log10(max(power,realmin))); xlim([0 50]); grid on;
    xlabel('Hz'); ylabel('Mean trial power (dB, source units)'); title('Per-trial spectra, then average');
    subplot(2,2,4);
    plot(time,mean(trials,1),'k'); grid on;
    xlabel('Time from face onset (s)'); ylabel('Source units');
    title(sprintf('Trial mean | Fp1 r=%.3f; Fp2 r=%.3f',r1,r2));
    print(fig,fullfile(root,sprintf('ICA%03d.png',ic)),'-dpng','-r140'); close(fig);
end
T=table((1:nc)',metrics(:,1),metrics(:,2),metrics(:,3),metrics(:,4), ...
    'VariableNames',{'component_one_based','source_variance','r_Fp1','r_Fp2','max_abs_frontal_correlation'});
writetable(T,fullfile(root,'component_metrics.csv'));
[~,order]=sort(T.max_abs_frontal_correlation,'descend');
writetable(T(order,:),fullfile(root,'frontal_correlation_ranking.csv'));
for start=1:16:nc
    fig=figure('Visible','off','Color','w','Position',[10 10 1200 1000]);
    for j=start:min(start+15,nc)
        subplot(4,4,j-start+1); topoplot(ICA_mm(:,j),h.chanlocs,'electrodes','on');
        title(sprintf('ICA%d',j));
    end
    print(fig,fullfile(root,sprintf('topographies_%03d_%03d.png',start,min(start+15,nc))),'-dpng','-r140'); close(fig);
end
fid=fopen(fullfile(root,'README.txt'),'w');
fprintf(fid,['All component numbers are native Letswave ONE-BASED.\n' ...
    'Fp1/Fp2 are scalp EEG proxies, NOT confirmed EOG. Correlation ranking is not a deletion list.\n' ...
    'Every component has a PNG; no component has been applied/excluded.\n' ...
    'ICA signs and numbering may change between independent fits; review this specific saved matrix.\n']); fclose(fid);
end
