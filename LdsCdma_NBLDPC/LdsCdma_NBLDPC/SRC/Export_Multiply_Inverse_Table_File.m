clear; clc;

qAry = 1024;

%% File Name
TableFileName =  sprintf('Arith.Table.GF.%d.txt', qAry);

A_Poly_ALL = primpoly(log2(qAry),  'all');
% A_Poly = A_Poly_ALL(1);
A_Poly = primpoly(log2(qAry));

global MULTIPLY_TABLE;
global INVERSE_TABLE;
global ADD_TABLE;

%% Generate the GF(q) table
NonBinary_LDPC_GenerateTable(qAry, A_Poly);

%% Export the Table to File
fid = fopen(TableFileName,'w');
fprintf(fid,'GF(%d) with Primitive Polynomial: %d. \n', qAry, A_Poly);
fprintf(fid,'Multiply Table:\n');
for j=1:qAry
    for i=1:qAry
        fprintf(fid,'%d ',MULTIPLY_TABLE(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'Add Table:\n');
for j=1:qAry
    for i=1:qAry
        fprintf(fid,'%d ',ADD_TABLE(i,j));
    end
    fprintf(fid,'\n');
end

fprintf(fid,'Inverse Table:\n');
fprintf(fid,'%d ',0);
for i=1:qAry-1
    fprintf(fid,'%d ',INVERSE_TABLE(i));
end
fprintf(fid,'\n');
fclose(fid);