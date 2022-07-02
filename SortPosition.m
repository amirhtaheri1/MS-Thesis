function pos = SortPosition(tmp)
% Sorting the position of the FIS variables used in parEPSO_OP.
% tmp is the position vector.
% pos is the sorted position vector approporiated to be used in the fis in
% detail.m.
pos = [sort(tmp(1:7)), sort(tmp(8:17)), ...
       sort(tmp(18:24)), sort(tmp(25:34)), sort(tmp(35:44)),...
       sort(tmp(45:51)), sort(tmp(52:61)), sort(tmp(62:71)),...
       sort(tmp(72:78))];
end

