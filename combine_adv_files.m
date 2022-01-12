function combine_adv_files(fn_list)
% fn_list: cell array containing filenames of the ADV .mat files to combine

for i = 1:length(fn_list)
    if i == 1
        load([fn_list{i} '.mat']);
    else
        partial_file = load([fn_list{i} '.mat']);
        Data.Profiles_TimeStamp = [Data.Profiles_TimeStamp; partial_file.Data.Profiles_TimeStamp];
        Data.Profiles_VelX = [Data.Profiles_VelX; partial_file.Data.Profiles_VelX];
        Data.Profiles_VelY = [Data.Profiles_VelY; partial_file.Data.Profiles_VelY];
        Data.Profiles_VelZ1 = [Data.Profiles_VelZ1; partial_file.Data.Profiles_VelZ1];
        Data.Profiles_VelZ2 = [Data.Profiles_VelZ2; partial_file.Data.Profiles_VelZ2];
        Data.Profiles_CorBeam1 = [Data.Profiles_CorBeam1; partial_file.Data.Profiles_CorBeam1];
        Data.Profiles_CorBeam2 = [Data.Profiles_CorBeam2; partial_file.Data.Profiles_CorBeam2];
        Data.Profiles_CorBeam3 = [Data.Profiles_CorBeam3; partial_file.Data.Profiles_CorBeam3];
        Data.Profiles_CorBeam4 = [Data.Profiles_CorBeam4; partial_file.Data.Profiles_CorBeam4];
        Data.Profiles_SNRBeam1 = [Data.Profiles_SNRBeam1; partial_file.Data.Profiles_SNRBeam1];
        Data.Profiles_SNRBeam2 = [Data.Profiles_SNRBeam2; partial_file.Data.Profiles_SNRBeam2];
        Data.Profiles_SNRBeam3 = [Data.Profiles_SNRBeam3; partial_file.Data.Profiles_SNRBeam3];
        Data.Profiles_SNRBeam4 = [Data.Profiles_SNRBeam4; partial_file.Data.Profiles_SNRBeam4];
        Data.Profiles_AmpBeam1 = [Data.Profiles_AmpBeam1; partial_file.Data.Profiles_AmpBeam1];
        Data.Profiles_AmpBeam2 = [Data.Profiles_AmpBeam2; partial_file.Data.Profiles_AmpBeam2];
        Data.Profiles_AmpBeam3 = [Data.Profiles_AmpBeam3; partial_file.Data.Profiles_AmpBeam3];
        Data.Profiles_AmpBeam4 = [Data.Profiles_AmpBeam4; partial_file.Data.Profiles_AmpBeam4];
    end
end

save([fn_list{1}(1:end-2) '.mat'],'Config','Data','-v7.3');