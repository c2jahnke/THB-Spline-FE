function errExport(fileName,H1err,L2err,Enerr,obj)
 % function for case study, write error to file
 fileID = fopen(fileName,'a');

 fmt = '\n %4f %4.16f %4.16f %4.16f %4.16f  %5f %4f';
 fprintf(fileID,fmt, [obj.levelBas{1}.N, H1err, L2err, Enerr,  obj.nOF obj.levelBas{1}.refArea]);
 disp("Error succesfully exported.");
 fclose(fileID);
end