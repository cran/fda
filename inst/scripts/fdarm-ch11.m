%
%
% Ramsay, Hooker & Graves (2009)
% Functional Data Analysis with R and Matlab (Springer)
%
% ch. 11  Functional Models and Dynamics
%

%  Set up some strings for constructing paths to folders.
%  These strings should be modified so as to provided access
%  to the specified folders on your computer.

%  Path to the folder containing the Matlab functional data analysis
%  software

fdaMPath = 'c:/Program Files/MATLAB/R2009a/fdaM';

%  Path to the folder containing the examples

examplesPath = [fdaMPath,'/examples'];

addpath(examplesPath)

%
% Section 11.1  Introduction to Dynamics
%
%  Figure 11.1

%
% Section 11.2 Principal Differential Analysis for Linear Dynamics
%
%  (no computations in this section)

%
% Section 11.3 Principal Differential Analysis of the Lip Data
%
% Figure 11.2

bwtCell = cell(2,1);
bwtCell{1} = fdPar(lipbasis,2,0);
bwtCell{2} = fdPar(lipbasis,2,0);

pdaStr = pda_fd(lipfd, bwtCell);

plot_pda(pdaStr)

dfd = 0.25*pdaStr.bwtCell[[2]]$fd^2 ...
         - pdaStr.bwtCell[[1]]$fd
             
dfd$fdnames = Cell(’time’,’rep’,’discriminant’)

% Figure 11.3 ????



% Figure 11.4
pda_overlay(pdaList)



%
% Section 11.4 PDA of the Handwriting Data
%

xfdCell = cell(3,1);
xfdCell{1} = fdafd(:,1);
xfdCell{2} = fdafd(:,2);
xfdCell{3} = fdafd(:,3);

pdaPar = fdPar(fdabasis,2,0);

pdaParCell = cell(2,1);
pdaParCell{1} = pdaPar;
pcaParCell{2} = pdaPar;

bwtCell = cell(2,2);
for i=1:2
    for j=1:2
        bwtCell{i,j} = pdaParCell;
    end
end
          
pdaCell = pda_fd(xfdCell, bwtCell);


eigen_pda(pdaCell)

% Figure 11.5









%
% Section 11.5 Registration and PDA
%
lipregCell = landmarkreg(lipfd, as.matrix(lipmarks),
    lipmeanmarks, WfdPar)
Dlipregfd = register.newfd(deriv_fd(lipfd,1),
    lipregCell$warpfd)
D2lipregfd = register.newfd(deriv_fd(lipfd,2),
    lipregCell$warpfd)
xfdCell = Cell(-Dlipregfd,-lipregCell$regfd)
lipregpda = fRegress( D2lipregfd, xfdCell, bwtCell)



%
% Section 11.6 Details for pda_fd, eigen_fd, pda_overlay
%              and register_newfd
%



%
% Section 11.7 Some Things to Try
%
% (exercises for the reader)

%
% Section 11.8  More to Read
%
