colnames=[{'fn'},fieldnames(measurement_data)'];
ds=[subjectlist',squeeze(struct2cell(measurement_data))'];
ds_f=vertcat(colnames,ds);
fid = fopen('/Users/jasongallant/Desktop/PKINGS_EODS/output/measurement_data.csv','wt');

if fid>0
 for k=1:size(C,1)
     if k==1
        fprintf(fid,'%s,%s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',ds_f{k,:});
     else
        fprintf(fid,'%s,%f,%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',ds_f{k,:});
     end
 end
  fclose(fid);
end
