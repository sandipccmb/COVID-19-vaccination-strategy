clear all; 

[num, txt, raw] = xlsread('Age_data_COVID_India2.xlsx');

[r,c] = find(strcmp(raw,'N'));
mat1 = cell2mat(raw(r+1:r+9,c+1:c+9)');

[r,c] = find(strcmp(raw,'Infected'));
mat2 = cell2mat(raw(r+1:r+9,c+1:c+9)');

[r,c] = find(strcmp(raw,'Demographics'));
demo = cell2mat(raw(r+1,c+1:c+9));
ppop = demo/sum(demo);

% --- Aggregate down the rows
mul = zeros(3,size(mat1,1));

mul(1,[1:3]) = 1;   % 0-29 y
mul(2,[4:6]) = 1;   % 30- 64 y
mul(3,[7:end]) = 1; % >64                                                % This matrix specifies which age categories to club together - tailor to your needs

% Used in the 1st submission of BMJ Open
% mul(1,[1:3]) = 1;   % 0-29 y 
% mul(2,[4:5]) = 1;   % 30-49 y
% mul(3,[6:end]) = 1; % >50                                                % This matrix specifies which age categories to club together - tailor to your needs

mat1a = mul*mat1;                                                          % Matrix based on reported contacts
mat2a = mul*mat2;                                                          % Matrix based on infections amongst contacts

% --- Aggregate across the columns
mat1b = (mul*(mat1a.*repmat(ppop,3,1))')';
mat2b = (mul*(mat2a.*repmat(ppop,3,1))')';

% Check that the two matrices are roughly in proportion
props = sum(mat2b,1)./sum(mat1b,1);

% --- Use the contact matrix, normalising so that mean value is 1
contact_mat = mat1b/mean(mat1b(:))

