%% Tutorial 04: Time-implicit multi-phase level set segmentation
%  Martin Rajchl, Imperial College London, 2015
%
%   [1] Rajchl, M.; Baxter, JSH.; Bae, E.; Tai, X-C.; Fenster, A.; 
%       Peters, TM.; Yuan, J.;
%       Variational Time-Implicit Multiphase Level-Sets: A Fast Convex 
%       Optimization-Based Solution
%       EMMCVPR, 2015.


clear all; close all;

% include max-flow solver
addpath(['..', filesep, 'maxflow']);
addpath(['..', filesep, 'lib']);

% 1. Load an image to segment
load('../data/medical_imgs.mat','brain_axial');
img = brain_axial;

% 2. Normalize the image intensity to [0,1]:
img = single(img);
img_n = (img - min(img(:))) / (max(img(:)) - min(img(:)));

[sx, sy] = size(img_n);
numberOfLabels = 3;
labelIDs = 1:numberOfLabels;

% 3. Create [numberOfLabels] initial regions
regions = zeros(size(img_n),'like', img_n); 

sx2 = idivide(sx,uint8(2));
sy2 = idivide(sy,uint8(2));

regions(1:sx2,1:sy2) = 1;
regions(sx2:end,1:sy2) = 2;
regions(1:end,sy2:end) = 3;

% 4. Construct a Potts graph according to [1]:
Ct = zeros(sx,sy,numberOfLabels);

% allocate alpha(x), the regularization weight at each node (x,label)
alpha = zeros(sx,sy,numberOfLabels);

% 5. Set up parameters and start level set iterations:
maxLevelSetIterations = 10; % number of maximum time steps
tau = 250; % speed parameter
w1 = 0.9; % weight parameter for intensity data term
w2 = 0.1; % weight parameter for the speed data term

for t=1:maxLevelSetIterations
        
    % 6. Compute a data term for each label according to [1]
    for i=1:numberOfLabels
        
        l = labelIDs(i);
        currRegion = regions == l;
        
        % 7. Compute a speed data term based on the sub-region extent as [1]:
        d_speed_in = -bwdist(1-currRegion,'euclidean').*currRegion;
        d_speed_out = bwdist(currRegion,'euclidean').*(1 - currRegion);
        
        d_speed = (d_speed_in + d_speed_out)./tau;
                
        % 8. Compute an simple example intensity data term based on the 
        % L1 distance to the mean of the current region
        m_int_inside = mean(mean(img_n(currRegion == 1)));
               
        d_int_inside = abs(img_n - m_int_inside);
                
        d_int = d_int_inside;
        
        % 9. Weight the contribution of both costs and assign them as sink 
        % capacities Ct in the graph. Note, that in the Potts 
        % configuration the source flows are unconstrained.
        Ct(:,:,i) = w1.*d_int + w2.*d_speed;
        
    end
    
    % 10. Assign a regularization weight (equivalent to pairwise terms) for
    % each node (x, label). Here we employ a constant regularization weight 
    % alpha. The higher alpha is, the more smoothness penalty is assigned.
    alpha = 0.05.*ones(sx,sy,numberOfLabels);
    
    % 11. Set up the parameters for the max flow optimizer:
    % [1] graph dimension 1
    % [2] graph dimension 2
    % [3] number of labels
    % [4] number of maximum iterations for the optimizer (default 200)
    % [5] an error bound at which we consider the solver converged (default
    %     1e-5)
    % [6] c parameter of the multiplier (default 0.2)
    % [7] step size for the gradient descent step when calulating the spatial
    %     flows p(x) (default 0.16)
    pars = [sx; sy; numberOfLabels; 200; 1e-5; 0.75; 0.16];
    
    % 12. Call the max flow optimizer with Ct, alpha and pars to 
    % obtain the continuous labelling function u, the convergence over 
    % iterations (conv), the number of iterations (numIt) and the run 
    % time (time):
    [u, conv, numIt, time] = asetsPotts2D(Ct, alpha, pars);
    
    % 13. Majority vote the continuous labelling function u to obtain a 
    % discrete segmentation result
    [uu,regions] = max(u, [], 3);
    
    % 14. Reassign the existing label IDs and recalulate the number of
    % labels
    labelIDs = unique(regions);
    numberOfLabels = length(labelIDs);
    
    % visualize the data terms Ct, contours and discretized labeling
    % functions:
    close all;
    figure();
    col = 'rgbcmyk';
    maxCt = max(Ct(:));
    minCt = min(Ct(:));
    for i=1:numberOfLabels
        subplot(3,numberOfLabels,i); imshow(Ct(:,:,i),[minCt maxCt]); title(['Ct(',num2str(i),')']);
        subplot(3,numberOfLabels,i+numberOfLabels); imshow(img_n,[]); 
        hold on; contour(regions == i, col(i)); title(['R_',num2str(i)]);
        subplot(3,numberOfLabels,i+2*numberOfLabels); imshow(regions == i,[]); title(['u_',num2str(i)]);
    end
    drawnow();
        
end



