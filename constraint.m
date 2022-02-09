function [Cineq,Ceq] = constraint(x)
global  colthickness slength ...
    tlthickness blthickness fid fid1 fid4 fid5
slength=0.8;
colthickness=0.1;
tlthickness=0.01;
blthickness=0.01;
% Set the constraints.
% Maximum absolute value of 
% vertical displacement and crushing force.
Dmaxver=40;
Fpeak=30000;
% Write the parameterized input file from Matlab.
fileID=fopen('honeycomb_sandwich.txt','w');
fprintf(fileID,'%s\n',...
    [num2str(slength),' ',num2str(x(:,1)),' ',...
    num2str(x(:,2)),' ',num2str(colthickness),' ',....
    num2str(tlthickness),' ',num2str(blthickness)]);
fclose(fileID);
fprintf('%d\n',x(:,1));
fprintf('%d\n',x(:,2));
% Run the Python script from Matlab through Abaqus.
system(['abaqus cae noGUI=Crush_sandwich.py']);
pause(10)
while exist('HoneyCombStructureJob1.lck','file')==2
    pause(0.1)
end
% Write the maximum nodal displacement.
fID=fopen('honeycombstructure_output_1.txt','r');
formatSpec='%f';
OutData=fscanf(fID,formatSpec);
fclose(fID)
% Write the maximum crushing force.
maxNodDisplY1=OutData(1,1)
fID=fopen('honeycombstructure_output_2.txt','r');
formatSpec='%f';
OutData=fscanf(fID,formatSpec);
fclose(fID)
Force=OutData(1,1)
% Constraint inequality matrix.
Cineq = [maxNodDisplY1-Dmaxver,Force-Fpeak];
Ceq = [];
fprintf(fid,'%4.2f\n',Force);
fprintf(fid1,'%f\n',maxNodDisplY1);
fprintf(fid4,'%f\n',(x(:,1)));
fprintf(fid5,'%4.2f\n',(x(:,2)));
end