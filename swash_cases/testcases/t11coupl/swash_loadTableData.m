function data = swash_loadTableData(file,np,nc)
    fid = fopen(file);

    if ~(fid>0)
        display(['Could not open file: ',file]);
        data = cell(0,0);
        return
    end

    names = {'time'};

    if nc>1
        for ic=2:nc
            names{ic} = ['Col',num2str(ic,'%02i')];
        end
    end

    data = cell(nc+1,1);
    data{1} = names;

    str = '';
    for ic=1:nc
        str = [str,'%f'];
    end
    res = cell2mat(textscan(fid,str));

    N = length(res(:,1));

    NT = floor(N/np);

    time = reshape(res(1:(NT*np),1),np,NT);
    data{2} = time(1,:)';
    for ic=2:nc
        dat = reshape(res(1:(NT*np),ic),np,NT);
        data{ic+1} = dat';
    end

    fclose(fid);
