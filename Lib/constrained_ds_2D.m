function xd = constrained_ds_2D(dsHandle, x, xGoal, constrCoef, constrMode)
    if(nargin<5)
        constrMode = 'zcbf'; % Zeroing control barrier function
    end

    % Original DS velocity
    xd = dsHandle(x - xGoal);
    % Compute constraint blouds
    numC = size(constrCoef,1);
    % Dimension of the state-space
    spaceDim = size(xd,1);
    
    unom = zeros(spaceDim, 1);
    
    %% Constrain DS dynamics
    if(strcmp(constrMode,'zcbf')) % Zeroing control barrier function
        A = [];
        b = [];
        gamma = 10;
        for i=1:numC
            n = constrCoef(i,1:spaceDim)';
            h = n'*x+constrCoef(i,spaceDim+1);

            if(n'*xd + gamma*h <= 0.0)
                A = [A; -n'];
                b = [b; n'*xd + gamma*h];
            end
        end

        if((~isempty(A)) && (rank(A)==size(A,1)))
            unom = pinv(A)*(b - A*unom) + unom;
        end
    elseif(strcmp(constrMode,'rcbf')) % Reciprocal control barrier function
        A = [];
        b = [];
        gamma = 100;  

        for i=1:numC
            n = constrCoef(i,1:spaceDim)';
            h = n'*x+constrCoef(i,spaceDim+1);
            
            % B = -ln(h/(1+h)) is a valid rcbf
            LfB = -n'*xd/(h*(1+h)); % Lie derivative of the barrier wrt f(x)
            LgB = -n'/(h*(1+h));

%             if(isnan(LfB))
%                 LfB = Inf;
%             end
            if(LfB - gamma*h >= 0.0)
                A = [A; LgB];
                b = [b; (gamma*h-LfB)];
            end
        end

        if((~isempty(A)) && (rank(A)==size(A,1)))
              %options = optimoptions('quadprog','Display','off');
%               unom = quadprog(eye(2),unom, A, b,[],[],[],[],[],options);
            
            unom = pinv(A)*(b - A*unom) + unom;
        end
    elseif(strcmp(constrMode,'ic'))   
        %% Invariance Control 
        %  Wolff and Buss "Invariance control design for constrained
        %  nonlinear systems", IFAC, 2005
        
        Phi = zeros(1,numC);
        ui  = zeros(1,numC);
        Ip = [];
        Im = [];
        % IC assumes negative barriers (i.e. h(x)<=0 is admissible)
        constrCoef = -constrCoef;  
        
        %constrCoef(end,3) = -constrCoef(end,3);
        for i=1:numC
            % Only for linear constraints with relatve degree = 1 (gamma = 0)
            Phi(i) = constrCoef(i,1)*x(1)+constrCoef(i,2)*x(2)+constrCoef(i,3);
        end
        PhiAll = max(Phi);
        if(PhiAll > -0.01)
            for i=1:numC
                % Compute a controller for each Phi
                n = constrCoef(i,1:2)';
                LfH = n'*xd;
                LgH = sum(n);
                ui(i) = -LfH / LgH;

                if(LgH >= 0)
                    Ip = [Ip i];
                else
                    Im = [Im i];
                end
            end 
            up = min([ui(Ip),Inf]);
            um = max([ui(Im),-Inf]);
            
            if(unom(1) > up)
                unom(1) = up;
            elseif(unom(1) < um)
                unom(1) = um;
            end
            
            if(unom(2) > up)
                unom(2) = up;
            elseif(unom(2) < um)
                unom(2) = um;
            end
        end
%         A = [];
%         b = [];
%         Phi = zeros(1,numC);
%         % IC assumes negative barriers (i.e. h(x)<=0 is admissible)
%         constrCoef = -constrCoef;  
%         for i=1:numC
%             % Only for linear constraints with relatve degree = 1 (gamma = 0)
%             Phi(i) = constrCoef(i,1)*x(1)+constrCoef(i,2)*x(2)+constrCoef(i,3);
%             if(Phi(i)> -0.01)   
%             	n = constrCoef(i,1:2)';
%                 
%                 A = [A; n'];
%                 b = [b; -n'*xd];
%             end
%         end
%         
%         if((~isempty(A)) && (rank(A)==size(A,1)))
%             % Works
%             options = optimoptions('quadprog','Display','off');
%             unom = quadprog(eye(2),unom, A, b,[],[],[],[],[],options);
%             
%             %unom = pinv(A)*(A*unom-b) + unom;
%         end
    end

    % Compute Invariance Control
%     A = diag([1 1 -1 -1]);
%     b =  [xd(1,1), xd(2,1), -xd(1,1), -xd(2,1)]';
%     b = - b;
% 
%     uds = A\b;
%     uM = [max(uds(1),-inf); max(uds(2),-inf)];
%     um = [min(uds(3),inf); min(uds(4),inf)];

%%  Working solution   
%     Phi = zeros(1,numC);
%     yd = Phi;
%     for i=1:numC
        % Only for linear constraints with relatve degree = 1 (gamma = 0)
%         Phi(i) = constrCoef(i,1)*x(1)+constrCoef(i,2)*x(2)+constrCoef(i,3);
%         yd(i)  = constrCoef(i,2)*xd(2)+constrCoef(i,1)*xd(1);
%         if(Phi(i)>=-0.01 && yd(i)>=0)
%             if(constrCoef(i,1)==0)
%                 unom(2) = -abs(constrCoef(i,2))*xd(2);
%             elseif(constrCoef(i,2)==0)
%                 unom(1) = -abs(constrCoef(i,1))*xd(1);
%             else
%                 unom(1) = -min(1,exp(500*Phi(i)))*(constrCoef(i,1)*xd(1)+constrCoef(i,2)*xd(2))/constrCoef(i,1);
%                 unom(2) = -min(1,exp(500*Phi(i)))*(constrCoef(i,1)*xd(1)+constrCoef(i,2)*xd(2))/constrCoef(i,2);
%                 
%                 unom(1) = -constrCoef(i,1)*xd(2); %-(constrCoef(i,1)*xd(1)+constrCoef(i,2)*xd(2))/constrCoef(i,1);
%                 unom(2) = constrCoef(i,2)*xd(1); %-(constrCoef(i,1)*xd(1)+constrCoef(i,2)*xd(2))/constrCoef(i,2);
%             end
%         end
%     end
%     xd = xd + unom;

    %% Test 1 - Works
%     M = eye(2);
%     Phi = zeros(1,numC);
%     yd = Phi;
%     for i=1:numC
        % Only for linear constraints with relatve degree = 1 (gamma = 0)
%         Phi(i) = constrCoef(i,1)*x(1)+constrCoef(i,2)*x(2)+constrCoef(i,3);
%         yd(i)  = constrCoef(i,2)*xd(2)+constrCoef(i,1)*xd(1);
%         if(Phi(i)>=-0.01 && yd(i)>=0)
%             if(constrCoef(i,1)==0)
%                 M(2,2) = 1-abs(constrCoef(i,2));
%             elseif(constrCoef(i,2)==0)
%                 M(1,1) = 1-abs(constrCoef(i,1));
%             else
%                 x0 = x(1);
%                 y0 = -(constrCoef(i,1)*x(1)+constrCoef(i,3))/constrCoef(i,2);
%                 
%                 x1 = 0;
%                 y1 = -constrCoef(i,3)/constrCoef(i,2);
%                 
%                 t = 1;
%                 if(abs(xd(1))>1e-8)
%                     M(1,1) = t*(x1-x0)/xd(1);
%                 end
%                 if(abs(xd(2))>1e-8)
%                     M(2,2) = t*(y1-y0)/xd(2);
%                 end
%             end
%         end
%     end
%     xd = M*xd;

    %% Test 2 - Works
%     Phi = zeros(1,numC);
%     yd = Phi;
%     for i=1:numC
        % Only for linear constraints with relatve degree = 1 (gamma = 0)
%         Phi(i) = constrCoef(i,1)*x(1)+constrCoef(i,2)*x(2)+constrCoef(i,3);
%         yd(i)  = constrCoef(i,2)*xd(2)+constrCoef(i,1)*xd(1);
%         if(Phi(i)>=0.0 && yd(i)>=0.01)     
%             t = 1;
%             if(constrCoef(i,1)==0)
%                 unom(2) = -xd(2)/abs(constrCoef(i,2));
%             elseif(constrCoef(i,2)==0)
%                 unom(1) = -xd(1)/abs(constrCoef(i,1));
%             else                
%                 x0 = x(1);
%                 y0 = -(constrCoef(i,1)*x(1)+constrCoef(i,3))/constrCoef(i,2);
%                 
%                 x1 = 0;
%                 y1 = -constrCoef(i,3)/constrCoef(i,2);
%                 
%                 u = [x1-x0; y1-y0];
%                 u = u / norm(u);
%                 
% %                 A = [constrCoef(i,1) constrCoef(i,2)];
% %                 b = -[constrCoef(i,1) constrCoef(i,2)]*xd;
% %                 
% %                 unom = pinv(A)*b;
% %                 
% %                 tmp_ = unom(2);
% %                 unom(2) = unom(1);
% %                 unom(1) = -unom(2);
%             
%                 unom(1) = -xd(1) + t*u(1); %-constrCoef(i,3)-xd(1)/abs(constrCoef(i,1)); %
%                 unom(2) = -xd(2) + t*u(2); %-constrCoef(i,3)-xd(2)/abs(constrCoef(i,2));
%             end
%         end
%     end



    xd = xd + unom;
end