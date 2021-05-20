function add2path
paths = regexp(genpath(cd),';','split');
nongit = cellfun(@regexp,paths,repmat({'\.'},[1 size(paths,2)]),repmat({'match'},[1 size(paths,2)]),'UniformOutput',false);
paths = paths(cellfun(@isempty,nongit));
addpath(paths{:});
savepath ;
end