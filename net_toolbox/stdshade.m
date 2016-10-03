function stdshade(amatrix,alpha,acolor,typem,smth,F)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23

if ~exist('acolor','var')
    acolor='r'; 
end

if ~exist('F','var')
    F=1:size(amatrix,2);
end
if ne(size(F,1),1)
    F=F';
end

if ~exist('smth','var');
    smth=1;
end  

if ~exist('typem','var')
    typem='mean'; 
end

switch typem
    case 'mean'
        amean=smooth(nanmean(amatrix),smth)';
    case 'median'
        amean=smooth(nanmedian(amatrix),smth)';
end

astd=nanstd(amatrix); % to get std shading
astd=nanstd(amatrix)/sqrt(size(amatrix,1)); % to get sem shading

if exist('alpha','var')==0 || isempty(alpha) 
    fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
%     acolor='k';
else fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');    
end

if ishold==0
    check=true; else check=false;
end

hold on;plot(F,amean,'color',acolor,'linewidth',.75); %% 'k' was acolor; linewidht was 1.5; change color or linewidth to adjust mean line

if check
    hold off;
end

end



