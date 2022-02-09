function F = objective(x)
persistent n_calls;%<-this will be incrementing each time the function is called
if isempty(n_calls)
    n_calls=0;
end
global variable_tracking_array fid2 fid3%<-this is so you can see the global variable
fID=fopen('honeycombstructure_output.txt','r');
formatSpec='%f';
OutData=fscanf(fID,formatSpec);
fclose(fID)
% Calculate the maximum vertical nodal displacements
E=OutData(1,1)
f1=-E



fID=fopen('honeycombstructure_output_3.txt','r');
formatSpec='%f';
OutData=fscanf(fID,formatSpec);
fclose(fID)
% Calculate the maximum vertical nodal displacements
Weight=OutData(1,1)
f2=Weight-2000;
fprintf('%d\n',f1)
fprintf('%d\n',f2)

F = [f1,f2];

 n_calls=n_calls+1;%<-increments each time the objective_function completes
 variable_tracking_array(n_calls,1)=f1;%<-you can supply any local variable name, not just variable_1 :-)
 variable_tracking_array(n_calls,2)=f2;
fprintf(fid2,'%4.2f\n',E);
fprintf(fid3,'%4.2f\n',f2);

end