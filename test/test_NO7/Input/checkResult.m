function success = checkResult
% automatic testscript for test case 7

success = true;

try

    % check if xml file was loade and interpreted correctly
    if ~exist(fullfile('tmp','PO320mg','applicationProtocol.mat'),'file')
        success = false;
        disp('failed in loading xml simulation file')
        return
    end
    
    load(fullfile('tmp','PO320mg','applicationProtocol.mat'));
    if isempty(ApplicationProtocol)
        success = false;
        disp('failed in loading xml simulation file')
        return
    end
    
    % check nonmem file read in
    if ~exist(fullfile('tmp','PO320mg','dataTp.mat'),'file')
        disp('failed in loading data file')
        return
    end
    
    load(fullfile('tmp','PO320mg','dataTp.mat'));
    
    % check for number of individuals
    nInd = size(DataTP,2);
    if nInd~=12
        disp(sprintf('nonmem load: 12 individuals expected, but there were %d',nInd));
        success = false;
        return
    end

    % check for fieldnames Data TP
    fnNorm = {'stud','sid','age','bmi','dose', 'gender', 'hght', 'wght'};
    fn = fieldnames(DataTP);
    if ~all(ismember(fn,fnNorm)) || ~all(ismember(fnNorm,fn))
        
        disp('nonmem load: fieldnames of DataTP ar unexpected');
        success = false;
        return
    end

    % check for number of individuals
    nInd = size(TP{1},2);
    if nInd~=12
        disp(sprintf('nonmem load: 12 individuals expected, but there were %d',nInd));
        success = false;
        return
    end

    
    % check for fieldnames TP
    fnNorm = {'time','dv','isLloq','lloq'};
    fn = fieldnames(TP{1});
    if ~all(ismember(fn,fnNorm)) || ~all(ismember(fnNorm,fn))
        
        disp('nonmem load: fieldnames of TP ar unexpected');
        success = false;
        return
    end
    
    % check values
    maxDV = max(cellfun(@(x) max(x), {TP{1}.dv}));
    if maxDV ~=  11.4000
        
        disp('nonmem load: maximal dv is unexpected');
        success = false;
        return
    end
    
    
catch exception
    
    disp(exception.message)
    success = false;
end

disp('automatic test was successful');

return
    