% script to generate numerical studies with 
% PoissGqThbMlGauss and PoissGqThbMlSin
 clc, clear all, close all
 % set basis file name bFname
 bFname = '/home/laptop/Documents/_mat_files/output/GaussHBP1P1NoDegRefFOO'
 Unif = false;
 % desP = [1,2,3,4,5];
 for p = 3:4
     
 fileName = strcat(bFname,num2str(p),'.txt') % add degree to file name
 if exist(fileName) == 2
     delete(fileName);
 end
 fileID = fopen(fileName,'a');
 fprintf(fileID,'N, H1err, L2err, Enerr,  obj.nOF refArea \n')
 fclose(fileID)
N = 0;
for k = 4 : 8
    
    N = 1*2^k
    PoissGqThbMlGauss(N,p,fileName);
end
if(Unif)
for k = 1 : 2
    N = N+64
     PoissGqThbMlGauss(N,p,fileName);
end
end
end