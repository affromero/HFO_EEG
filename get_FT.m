root = 'Result_Task_all';
Data_File = dir(root);
Data_File = {Data_File(3:end).name};
Data_mat  ={};
pacient   ={};

%Getting all file directions with data in .mat
for i=Data_File
    folder = dir(fullfile(root, cell2mat(i)));
    folder = {folder(3:end).name};
    for j=folder
        Data_mat{end+1} = fullfile(root, cell2mat(i), cell2mat(j));
        pacient{end+1}  = cell2mat(i);
    end
end

outD = 'Results_FT';
mkdir(outD);
for i = 1:numel(Data_mat)
    Result                = [];
    pacient_number        = Data_mat{i}(numel(root)+2+numel(pacient{i})+1:end-4);
    result_pacient_folder = fullfile(outD, pacient{i});
    mkdir(result_pacient_folder);
    outName               = ['HFO_' pacient_number '.mat'];
    outName               = fullfile(result_pacient_folder , outName);
    if exist(outName, 'file')
        fprintf('File %s already exist\n', outName);
        continue; 
    end
    out_conv = {};
    
    load(Data_mat{i})
    channels = fieldnames(data);

    for ch = 1:numel(channels)
        ch1 = char(channels(ch));
        fprintf('Detecting HFO over %s - %s channel\n', pacient_number, ch1);   
        if numel(data.(ch1))>(1024*240)
            data.(ch1) = data.(ch1)(end-(1024*240)+1:end);
        end
        [Result.(ch1), ~] = find_HFO_peaks(data.(ch1));
    end
    
    eval('Result');
    if isempty(Result), continue; end
    fprintf('...Saving Result in %s file\n\n', outName);
    save(outName, '-v7.3', 'Result');
end
    
    