function [Pfit, Pweigh, y_fit, X, Y] = MarijnCalcCorrelationStr(M)

OutDeg 	= sum(M > 0, 1)';
InDeg  	= sum(M > 0, 2);

Y       = OutDeg;
[X, Ix] = sort(InDeg);
Y       = Y(Ix);
[p,err] = polyfit(X,Y,1);
Pfit 	= p(1);
y_fit   = polyval(p,X,err); 

midp    = mean([mean(InDeg) mean(OutDeg)]);
PC      = sum(InDeg > midp & OutDeg > midp) + sum(InDeg < midp & OutDeg < midp);
AC      = sum(InDeg > midp & OutDeg < midp) + sum(InDeg < midp & OutDeg > midp);
Pweigh  = PC / AC - 1;