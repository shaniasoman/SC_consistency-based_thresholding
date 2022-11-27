
clear all;

%subject list
subjects = {'sub0086','sub0134','sub0152','sub0334','sub0412','sub0550','sub0575','sub0613','sub0616','sub0647','sub1110','sub1170','sub1264','sub1274','sub1714','sub1908','sub1938','sub2004','sub2031','sub2093','sub2105','sub2144','sub2498','sub2582','sub2654','sub2704','sub2819','sub4030','sub4264','sub4405','sub4545','sub4656','sub4789','sub4925','sub4970','sub5051','sub5180','sub5384','sub5459','sub5506','sub5515','sub5519','sub5579','sub5588','sub5662','sub5729','sub5856','sub5927','sub5932','sub5946','sub5965','sub6026','sub6047','sub6207','sub6279','sub6490','sub6569','sub6589','sub6618','sub6700','sub6751','sub6781','sub6812','sub6820','sub6909','sub6919','sub7019','sub7020','sub7074','sub7187'};
%concatenating structural conenctivity matrix of each subject
Ws=cat(3,sub0086,sub0134,sub0152,sub0334,sub0412,sub0550,sub0575,sub0613,sub0616,sub0647,sub1110,sub1170,sub1264,sub1274,sub1714,sub1908,sub1938,sub2004,sub2031,sub2093,sub2105,sub2144,sub2498,sub2582,sub2654,sub2704,sub2819,sub4030,sub4264,sub4405,sub4545,sub4656,sub4789,sub4925,sub4970,sub5051,sub5180,sub5384,sub5459,sub5506,sub5515,sub5519,sub5579,sub5588,sub5662,sub5729,sub5856,sub5927,sub5932,sub5946,sub5965,sub6026,sub6047,sub6207,sub6279,sub6490,sub6569,sub6589,sub6618,sub6700,sub6751,sub6781,sub6812,sub6820,sub6909,sub6919,sub7019,sub7020,sub7074,sub7187);
   
% group mean connectivity matrix
Wmean=mean(Ws,3);
% coefficient of variation across the group
Wcv=std(Ws,0,3)./Wmean; 
% removing NaN
Wcv(isnan(Wcv))=0;
% p - proportion of weights to preserve
p=prctile(Wcv,75);

subjectdir='C:\Users\Shania Jijo\OneDrive - Deakin University\Documents\PhD_Study1\Study 3\SC_Wave1\SC_allsub';

for i=1:length(subjects)
    sub=subjects{i};
    conmat=load(fullfile(subjectdir,sub,'.mat')); %load each subjects structural connectivity matrix
    %defining variables
    sname=conmat;
    %convert each connectivity matrix to upper triangular matrix
    t_sname = triu(sname,1);
    t_sname=threshold_arbmeasure(t_sname,-Wcv,p);
end

function W = threshold_arbmeasure(W, M, p)
%THRESHOLD_ARBMEASURE     Threshold edges ranked by given arbitrary measure
%
%   W_thr = threshold_arbmeasure(W, M, p);
%
%   This function "thresholds" the connectivity matrix by preserving a
%   proportion p (0<p<1) of the edges with the largest values of the given
%   measure M. All other weights, and all weights on the main diagonal
%   (self-self connections) are set to 0. 
%
%   Inputs: W,      weighted or binary connectivity matrix
%           M,      matrix of values of the measure by which to rank each
%                   node, preserving only edges with the largest values
%           p,      proportion of weights to preserve
%                       range:  p=1 (all weights preserved) to
%                               p=0 (no weights removed)
%
%   Output: W_thr,  thresholded connectivity matrix
%
%   James Roberts, QIMR Berghofer, 2014
%   Based directly on BCT's threshold_proportional.m by:
%     Mika Rubinov, U Cambridge,
%     Roan LaPlante, Martinos Center, MGH
%    (BCT: brain-connectivity-toolbox.net)

n=size(W,1);                                %number of nodes
W(1:n+1:end)=0;                             %clear diagonal

if isequal(W,W.');                          %if symmetric matrix
    W=triu(W);                              %ensure symmetry is preserved
    ud=2;                                   %halve number of removed links
else
    ud=1;
end

ind=find(W);                                %find all links
E=sortrows([ind M(ind)], -2);               %sort by M
en=round((n^2-n)*p/ud);                     %number of links to be preserved

W(E(en+1:end,1))=0;                         %apply threshold

if ud==2                                    %if symmetric matrix
    W=W+W.';                                %reconstruct symmetry
end

