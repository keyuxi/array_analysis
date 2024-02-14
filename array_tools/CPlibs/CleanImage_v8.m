function [image, zeroPix] = CleanImage_v7(image,maxVal,cutOff)
% initialize parameters
iSize = size(image);
backgroundEstimate = mode(image(:));
d = [1,0; -1,0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1]; % neighbor pixels

% get crap centers
[row,col] = find(image>=maxVal);
zeroPix = sub2ind(iSize,row,col);
status = 1;

% start fun
while status~=0;
    % remove boundry cases and initialize newPix
    row(row==1) = 2; row(row==iSize(1)) = iSize(1) - 1;
    col(col==1) = 2; col(col==iSize(2)) = iSize(2) - 1;
    
    % get neighbor pixels
    newPix = [];
    for i = 1:length(d);
        newPix = vertcat(newPix,sub2ind(iSize,d(i,1)+row,d(i,2)+col));
    end
    
    % get unique set to zero
    newPix = unique(newPix);
    newPix = newPix(image(newPix)>=cutOff);
    newPix = newPix(find(ismember(newPix,zeroPix)==0));
    
    % set saturated pix to maxVal
    image(newPix) = maxVal;
    image(zeroPix) = backgroundEstimate;
    
    % get all saturated regions
    [row,col] = find(image>=maxVal);
    zeroPix = unique(vertcat(zeroPix,sub2ind(iSize,row,col)));
    status = length(newPix);
end

% get surrounding indicies to remove residual crap [conservative]
for i = 1:2
    [row,col] = ind2sub(iSize,zeroPix);
    row(row==1) = 2; row(row==iSize(1)) = iSize(1) - 1;
    col(col==1) = 2; col(col==iSize(2)) = iSize(2) - 1;
    newPix = [];
    for i = 1:length(d);
        newPix = vertcat(newPix,sub2ind(iSize,d(i,1)+row,d(i,2)+col));
    end

    % set to bg and store
    newPix = unique(newPix);
    image(newPix) = backgroundEstimate;
    zeroPix = unique(vertcat(zeroPix,newPix));
end
end
