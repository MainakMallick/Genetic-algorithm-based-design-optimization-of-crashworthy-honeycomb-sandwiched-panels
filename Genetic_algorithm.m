% Open individual .txt file to monitor the i/o of
% every asscosiated metrics to monitor the
% improvements in each iteration.
global fid fid1 fid2 fid3 fid4 fid5
fid = fopen('Peak Force.txt','w');
fid1 = fopen('Nodal Displacement.txt','w');
fid2 = fopen('Energy absorbed.txt','w');
fid3 = fopen('Weight.txt','w');
fid4 = fopen('Wall thickness.txt','w');
fid5 = fopen('No of core cells.txt','w');
% All the syntaxed equalities/inequalities were
% passed as a null matrix as all of that will be 
% gathered from FE analysis.
Aineq = [0,0];
bineq = 0;
% The upper bound and the lower bound of intended
% sample space for the variable is provided.
lb = [0.002,10.0];
ub = [0.006,14.0];
fun = @objective;
nlcon = @constraint;
npts = 60; % The default is 60
opts_ga.ParetoSetSize = npts;
[x_ga1,fval_ga1,~,gaoutput1] = gamultiobj...
    (fun,2,Aineq,bineq,[],[],lb,ub,nlcon,opts_ga);
opts_ga = optimoptions('gamultiobj','Display',...
'off','PlotFcn','gaplotpareto','PopulationSize',npts);
disp("Total Function Count: " + psoutput2.funccount);
fprintf('%d\n',f1)
fprintf('%d\n',f2)
global variable_tracking_array;%<-create a global for tracking data
variable_tracking_array = zeros(1000,2);%<-in this case we're tracking 3 variables
