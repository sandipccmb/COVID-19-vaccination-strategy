[num, txt, raw] = xlsread('COMORBID_with_population.xlsx');

% Columns: 1-3 htn (with intvls) / 4-6 dm (with intvls) / 7-9 dm_htn (with intvls) / 10 Total
mat = permute(reshape(num,[4,2,36,10]),[1,4,2,3]);                         % Dims: 1.Age 2.Column 3.urbrur 4.State
