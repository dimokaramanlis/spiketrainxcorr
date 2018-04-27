function [XCmat] = xcorrSpikeTrains(spiketrains,maxLag,normflag)
%XCORRSPIKETRAINS Summary of this function goes here
%   Detailed explanation goes here
%   Input:
%       spiketrains: Nbins x Ncells matrix for all cells
%       maxLag: number of time bins for xcorr calculation 
%       normflag: 0 for raw xcorr, 1 for normalized ('coeff')
%   Output:
%       XCmat: 2*maxLag+1 x Ncells x Ncells 3D array
%------------------------------------------------------------------------
% made based on 
% https://stats.stackexchange.com/questions/120513/cross-correlation-for-very-sparse-binary-data
% written by Dimos, 27.04.2018
%------------------------------------------------------------------------

Ncells=size(spiketrains,2);
lag = -maxLag:maxLag;

%make spike trains sparse to increase matrix multiplication speed
sparsetrains = sparse(spiketrains); 

XCmat = zeros(length(lag),Ncells,Ncells);
for i=1:length(lag)
    %shift matrices
    if lag(i)>=0
        Y1=sparsetrains(1+lag(i):end,:); Y2=sparsetrains(1:end-lag(i),:);
    else
        Y1=sparsetrains(1:end+lag(i),:); Y2=sparsetrains(1-lag(i):end,:);
    end
    XCmat(i,:,:)=Y1'*Y2; %core calculation
end

%normalize so that the autocorrelations at zero lag equal 1
%similar to MATLAB's 'coeff'
if normflag; XCmat=normalizeXcorr(XCmat, maxLag); end

end

function Xnorm=normalizeXcorr(XCmat, maxLag)
    Ncells=size(XCmat,3);
    allZeroCorrs=diag(squeeze(XCmat(maxLag+1,:,:)));
    allDivisors=reshape(allZeroCorrs*allZeroCorrs',1,Ncells,Ncells);
    Xnorm=bsxfun(@rdivide,XCmat,sqrt(allDivisors));
end
