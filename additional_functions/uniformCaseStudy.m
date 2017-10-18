
 fileName = '/home/laptop/Documents/_mat_files/output/foo.txt';
 if exist(fileName) == 2
     delete(fileName);
 end
 fileID = fopen(fileName,'a');
 fprintf(fileID,'N, H1err, L2err, Enerr,  obj.nOF refArea \n')
 fclose(fileID)
 
N = 0;
for k = 3 : 6
    N = 1*2^k
    PoissGqThbMlSin(N);
end
% for k = 1 : 5
%     N = N+16;
%     PoissGqThbMlSin(N)
% end