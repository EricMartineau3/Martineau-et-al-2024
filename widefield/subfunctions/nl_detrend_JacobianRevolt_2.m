function [y_fit, coefs,gof]= nl_detrend_JacobianRevolt_2(x, data, bounds,chanName) 
% Fits an exponential function, y = a*x^b, for each pixel or for the average pixel, excluding the
% timepoints inside the bounds input (stim + expected reponse window) 
% Inputs:
%   x = time values as vector
%   data = video, inputed as XY-by-t matrix
%   bounds = exclusion bounds. Values to ignore during fit.
%   chanName = (optional) channel name, to display in waitbar. 
%
% Processing time has been optimized by splitting the image in 256 pixel
% bins. Jacobian method is faster than pixel-by-pixel loop at low number of
% fits(pixels) but solving the equation with higher number of fits get
% exponentially slower. Sweet-spot for optimal performance improvement has
% been empirically determined to be around 500. 
%
% Written by Eric Martineau and Antoine Malescot, based on Exemple from https://www.mathworks.com/matlabcentral/answers/431697-make-curve-fitting-faster

%%
%Prep
nPixels = size(data,1);
binSize = 256;
nBin = nPixels/binSize;
nanIdx = isnan(data(1,:)); %identify nan frames (movement);
data = data(:,~nanIdx); %remove NaN frames;
xx = x(:,~nanIdx); %remove NaN timepoints, but keep original timescale for detrend function later
indx = [find(xx>=bounds(1),1), find(xx>=bounds(2),1)];
data = data(:,[1:indx(1),indx(2):end]); %remove points outside of bound
xx_mat = repmat(xx(:,[1:indx(1),indx(2):end]),binSize,1);%remove points outside of bound

%Fit average image to estimate starting parameters
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', [1 -100 1 -0.1], 'Upper',[2^16 0 2^16 0.1],'TolX',1e-20,'TolFun',1e-20);
opts.Display = 'Off';
% opts.StartPoint = [116.607174931279 -0.260820095436706 5036.74930225973 0.00129765068984856];
opts.StartPoint = [10000 -0.01 1000 0.001];
[fitresult,gof,~] = fit( xx_mat(1,:)', mean(data,1)', 'exp2',opts);
x_m = coeffvalues(fitresult); %coeff of mean pixel trend
x0 = repmat(x_m,binSize,1);
clear x_m

%Define options
tolfun = 1e-6;
StepTol = 1e-4;
optsNJ=optimoptions('lsqcurvefit','SpecifyObjectiveGradient',true,'Jacobian','on','display','off', ...
    'Algorithm','trust-region-reflective','tolfun',tolfun,'StepTolerance',StepTol);
lb = repmat([-2^16,-Inf,-2^16,-Inf],binSize,1); %lower bound of each parameter
ub = repmat([2^16,Inf,2^16,Inf],binSize,1); %up bound of each parameter

%Curve fitting by splitting image in 256 pixel bins 
y_fit = zeros(size(data,1),size(x,2));
coefs = zeros(size(data,1),size(x0,2));
f = waitbar(0,strcat("Exponential trend fitting, ",chanName));
for jj = 1:nBin 
    idx = [(jj-1)*binSize+1,jj*binSize];
    coef = x0;
    [yy_fit,~]=fun_power(coef,repmat(x,binSize,1)); %Apply a and b coefficients to each pixel to get trend
    y_fit(idx(1):idx(2),:) = yy_fit;
    coefs(idx(1):idx(2),:) = coef;
    waitbar(jj/nBin,f);
end
clear yyfit coef jj
close(f)

%%%%%%%%%%%%%%%%%%Define functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [f,J] = fun_power(x,xdata) % Equation is Y = a*exp(b*x)+c*exp(d*x), where x is the parameters, x(:,1) = a, x(:,2) = b etc.; and xdata are the pixel values across time
        %xdata is a xy-by-t matrix; x is a xy-by-nparameters fit;
        npoints = size(xdata,2); %t
        nFits = size(xdata,1); %pixels
        nParam = size(x,2); 
    
        f = x(:,1)*ones(1,npoints) .* exp(xdata .* (x(:,2)*ones(1,npoints))) + x(:,3)*ones(1,npoints) .* exp(xdata .* (x(:,4)*ones(1,npoints)));
        J = spalloc(npoints*nFits, nParam*nFits, npoints*nFits*nParam);
    
        for i=1:nFits       % loop over fits, calculate 3*n derivatives for each fit
            fp=f(i,:);     % f of points for this fit (n)
            xp=xdata(i,:); % points for this fit (n)
            p=x(i,:);      % params for this fit (m)        
            J(i:nFits:end,        i)= exp(xp.*p(2));    % df/da     (may be negative)
            J(i:nFits:end,   nFits + i)=p(1).*xp.* exp(p(2).*xp); % df/db   (non-negative)
            J(i:nFits:end,  2*nFits + i)=exp(xp.*p(4));% df/dc 
            J(i:nFits:end,  3*nFits + i)=p(3).*xp.* exp(p(4).*xp);% df/dd 
            %got to wolframalpha.com to compute derivatives
        end     
    end
    
    function stop = myoutput(x,optimvalues,state) %output function to 
        global history;
        stop = false;
        if isequal(state,'iter')
            history = [history; x];
        end
    end
end
