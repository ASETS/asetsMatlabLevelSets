%% Tutorial 03: Time-implicit coupled level set segmentation
%  Martin Rajchl, Imperial College London, 2015
%
%   [1] Rajchl, M.; Baxter, JSH.; Bae, E.; Tai, X-C.; Fenster, A.; 
%       Peters, TM.; Yuan, J.;
%       Variational Time-Implicit Multiphase Level-Sets: A Fast Convex 
%       Optimization-Based Solution
%       EMMCVPR, 2015.
%
%   [2] Ukwatta, E.; Yuan, J.; Rajchl, M.; Qiu, W.; Tessier, D; Fenster, A.
%       3D Carotid Multi-Region MRI Segmentation by Globally Optimal 
%       Evolution of Coupled Surfaces
%       IEEE Transactions on Medical Imaging, 2013
%
%   [3] Rajchl M., J. Yuan, E. Ukwatta, and T. Peters (2012).
%       Fast Interactive Multi-Region Cardiac Segmentation With 
%       Linearly Ordered Labels. 
%       ISBI, 2012. pp.1409–1412.


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
numberOfLabels = 4;

% 3. Create [numberOfLabels] initial coupled regions
regions = ones(size(img_n),'like', img_n); % background
for i=2:numberOfLabels
    d = 20*i;
    regions(d:sx-d,d:sy-d) = i; % linearly ordered region ids i
end

% visualize the initial regions
col = 'rgbcmyk';
figure(); title('Initial regions');
imshow(img_n,[]);
hold on; 
for i=1:numberOfLabels
   contour(regions == i,col(i));
end
hold off;
drawnow();

% 4. Construct a linearly ordered graph according to [3]:
Ct = zeros(sx,sy,numberOfLabels);

% allocate alpha(x), the regularization weight at each node (x,label)
alpha = zeros(sx,sy,numberOfLabels-1);

% 5. Set up parameters and start level set iterations:
maxLevelSetIterations = 15; % number of maximum time steps
tau = 50; % speed parameter
w1 = 0.5; % weight parameter for intensity data term
w2 = 0.5; % weight parameter for the speed data term

for t=1:maxLevelSetIterations
    
    % 6. Compute a data term for each label according to [1]
    for i=1:numberOfLabels
        
        currRegion = regions == i;
        
        % 7. Compute a speed data term based on the sub-region extent as [1]:
        d_speed = ((1-currRegion).*bwdist(currRegion,'Euclidean'))./tau;
                
        % 8. Compute an simple example intensity data term based on the 
        % L1 distance to the mean of the current region
        m_int_inside = mean(mean(img_n(currRegion == 1)));
        m_int_outside =  mean(mean(img_n(currRegion == 0)));
        
        d_int_inside = abs(img_n - m_int_inside);
        d_int_outside = abs(img_n - m_int_outside);
        
        % 9. Weight the contribution of both costs and assign them as sink 
        % capacities Ct in the graph. Note, that in the Ishikawa 
        % configuration the source flows are unconstrained.
        Ct(:,:,i) = w1.*(d_int_outside + d_int_inside) + w2.*(d_speed);
    end
    
    % 10. Assign a regularization weight (equivalent to pairwise terms) for
    % each node (x, label). Here we employ a constant regularization weight 
    % alpha. The higher alpha is, the more smoothness penalty is assigned.
    alpha = 0.1.*ones(sx,sy,numberOfLabels-1);
    
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
    
    % 12. Call the Ishikawa max flow optimizer with Ct, alpha and pars to 
    % obtain the continuous labelling function u, the convergence over 
    % iterations (conv), the number of iterations (numIt) and the run 
    % time (time):
    [u, conv, numIt, time] = asetsIshikawa2D(Ct, alpha, pars);
    
    % 13. Threshold the continuous labelling function u to obtain a discrete
    % segmentation result
    regions = ones(sx,sy);
    for i=1:numberOfLabels-1
        regions = regions + (u(:,:,i) > 0.5);
    end
    
    
    % visualize the data terms Ct, contours and discretized labeling
    % functions:
    close all;
    figure();
    for i=1:numberOfLabels
        subplot(3,numberOfLabels,i); imshow(Ct(:,:,i),[]); title(['Ct(',num2str(i),')']);
        subplot(3,numberOfLabels,i+numberOfLabels); imshow(img_n,[]); 
        hold on; contour(regions == i, col(i)); title(['R_',num2str(i)]);
        subplot(3,numberOfLabels,i+2*numberOfLabels); imshow(regions == i,[]); 
    end
    drawnow();
    pause(0.5);
        
end


