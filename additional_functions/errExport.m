function errExport(H1err,L2err,Enerr,obj)
 fileName = '/home/laptop/Documents/_mat_files/output/AdaptErrp3Prec9.txt';
 
 fileID = fopen(fileName,'a');

 fmt = '\n %4f %4.9f %4.9f %4.9f %4.9f  %5f %4f \n';
 fprintf(fileID,fmt, [obj.levelBas{1}.N, H1err, L2err, Enerr,  obj.nOF obj.levelBas{1}.refArea]);
 disp("Error succesfully exported.");
 fclose(fileID);
end