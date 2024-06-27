function plot_linear_constraints(constrCoef, xLim, color)
    x1 = xLim.min(1):0.1:xLim.max(1);
    x2 = xLim.min(2):0.1:xLim.max(2);
    
    for i=1:size(constrCoef,1)
        if(constrCoef(i,1)==0)
            X1 = x1;
            X2 = repmat(-constrCoef(i,3)/constrCoef(i,2),1,length(x1));
        elseif(constrCoef(i,2)==0)
            X2 = x2;
            X1 = repmat(-constrCoef(i,3)/constrCoef(i,1),1,length(x2));
        else
            X1 = x1;
            X2 = -(constrCoef(i,1)/constrCoef(i,2))*X1-(constrCoef(i,3)/constrCoef(i,2));
        end
        
        plot(X1, X2, 'Color', color)
    end   
end