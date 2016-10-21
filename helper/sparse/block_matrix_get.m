function block = block_matrix_insert(B,rows,columns,i,j)
% Z = block_matrix_insert(B,rows,columns,block)
% Gets a submatrix from a block matrix at block position (i,j)
%
% Input:
%       B: Block matrix (see block_matrix.m for a description).
%    rows: Number of elements in each block row of the block matrix.
% columns: Number of elements in each block row of the block matrix.
%     i,j: Block indices denoting where the submatrix should be placed.
% Output:
%   block: The requested submatrix
%
% Example:
% >> magic(5)
%
% ans =
%
%    17    24     1     8    15
%    23     5     7    14    16
%     4     6    13    20    22
%    10    12    19    21     3
%    11    18    25     2     9
% >> B = block_matrix_get(magic(5),[2 3],[2 3],2,2)
% 
% B = block_matrix_get(magic(5),[2 3],[2 3],2,2)
% 
% B =
% 
%     13    20    22
%     19    21     3
%     25     2     9

assert(i <= length(rows),          'Index out of bounds.');
assert(j <= length(columns),       'Index out of bounds.');

r = [0 rows];
c = [0 columns];
offset_i0 = sum(r(1:i));
offset_in = sum(r(1:i+1));
offset_j0 = sum(c(1:j));
offset_jn = sum(c(1:j+1));
block = B((1+offset_i0):offset_in,(1+offset_j0):offset_jn);
