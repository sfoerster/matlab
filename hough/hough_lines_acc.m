function [H, theta, rho] = hough_lines_acc(BW, varargin)
    % Compute Hough accumulator array for finding lines.
    %
    % BW: Binary (black and white) image containing edge pixels
    % RhoResolution (optional): Difference between successive rho values, in pixels
    % Theta (optional): Vector of theta values to use, in degrees
    %
    % Please see the Matlab documentation for hough():
    % http://www.mathworks.com/help/images/ref/hough.html
    % Your code should imitate the Matlab implementation.
    %
    % Pay close attention to the coordinate system specified in the assignment.
    % Note: Rows of H should correspond to values of rho, columns those of theta.

    %% Parse input arguments
    p = inputParser();
    addParameter(p, 'RhoResolution', 1);
    addParameter(p, 'Theta', linspace(-90, 89, 180));
    parse(p, varargin{:});

    rhoStep = p.Results.RhoResolution;
    theta = p.Results.Theta;

    %% TODO: Your code here

    %rhoResolution = 1;
    %thetaResolution = 0.5;

    [M, N] = size(BW);
    maxRho = sqrt((M-1).^2 + (N-1).^2);
    
    rho = -maxRho:rhoStep:maxRho;
    %theta = -90:thetaResolution:89.999;
    
    H = zeros(length(rho),length(theta));
    
    total=sum(BW(:));
    current=0;
    
    y=0;
    for m=M:-1:1 % go down vertically (0 at top)
        y=y+1;
        for n=1:N % go from left to right (0 at left)
            x=n;
            
            if BW(m,n) %true means edge
                current=current+1;
                if mod(current,100) == 0
                    disp(['[hough_lines_acc] Progress: ' num2str(current) '/' num2str(total)])
                end
               
                % an upsampled approximation of a continuous rho as a
                % function of theta
                
                trueRho = x.*cosd(theta) + y.*sind(theta);
                
                invChRho = (length(rho) - 1)/(rho(end) - rho(1));
                binRhoIdx = round(invChRho*(trueRho - rho(1)) + 1);

                for nnRho=1:length(rho)
                    htidx = binRhoIdx == nnRho;
                    H(nnRho,htidx) = H(nnRho,htidx)+1;
                end
                
                
            end
            
            
        end
    end

end
