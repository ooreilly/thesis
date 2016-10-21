function B = block_matrix_add(B,rows,columns,i,j,block)
% Z = block_matrix_add(B,rows,columns,block)
% Adds a submatrix into a block matrix at block position (i,j)
%
% Input:
%       B: Block matrix (see block_matrix.m for a description).
%    rows: Number of elements in each block row of the block matrix.
% columns: Number of elements in each block row of the block matrix.
%     i,j: Block indices denoting where the submatrix should be placed.
%   block: The submatrix to insert into the block matrix. 
%          The size of the submatrix must match the number given in rows(i)
%          and columns(j).
%
% Example:
% >> Z = block_matrix([2 3],[2 3])
% >> B = block_matrix_insert(Z,[2 3],[2 3],1,2,ones(2,3));
% >> B = block_matrix_add(Z,[2 3],[2 3],1,2,ones(2,3));
%>> full(B)
%
%ans =
%
%     0     0     2     2     2
%     0     0     2     2     2
%     0     0     0     0     0
%     0     0     0     0     0
%     0     0     0     0     0


A = block_matrix_get(B,rows,columns,i,j);
B = block_matrix_insert(B,rows,columns,i,j,A+block);



