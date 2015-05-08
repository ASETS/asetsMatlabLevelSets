function [cost] = computeLogLikelihoodCost(img, lbl, epsilon)
%% Compute the log likelihood cost from given intensity samples

% cast to float
img = single(img);

% define the number of histogram bins
nBins = 256;

minI = min(img(:));
maxI = max(img(:));

% normalize image to 8 bit
img_n = ((img - minI)/ (maxI - minI)).*255.0;

% compute histogram
[binCounts] = histc(img_n(lbl == 1),linspace(0,255,nBins));

% normalize to compute the probabilities
binCounts = binCounts./sum(binCounts(:));

% compute LL
P = binCounts( uint16(img_n/ (256/nBins)) + 1);
cost = -log10(P  + epsilon);

end
