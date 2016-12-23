clear;
clc;

%% Configuration
Q = 64;
P = log2(Q);
NB_MAT_REPR_File =  sprintf('Mat.Repr.GF.%d.txt', Q);

%% Generate Primitive Matrics Rrepresentation 
A_Poly_ALL = primpoly(P, 'all');
A_Poly = primpoly(P);

coef_a = bitget(A_Poly, 1:P);
A_Matric = [[zeros(P-1,1),diag(ones(P-1,1))];coef_a];

%% Generate All Order Nonzero Element & Write to File
fid = fopen(NB_MAT_REPR_File, 'w');
fprintf(fid,'GF(%d) with Primitive Polynomial: %d \n', Q, A_Poly);
alpha = gf(2, P, A_Poly); % since the poly of alpha is 00...010 = 2
for k = 0:Q-2
    a_k = gf_transform_2_int(alpha^k, P, A_Poly);
    A_K = mod(A_Matric^k,2);
    fprintf(fid,'A^%d --> order: %d\tpoly: %d\n', k, k, a_k);
    for i = 1:P
        for j = 1:P
            fprintf(fid,'%d ',A_K(i,j));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
