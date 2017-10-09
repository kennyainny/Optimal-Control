
    %movement from point to point
    %trajVars = {6, 2, 3};
    %moves(1).direction = {0, 0, 0}; %x, y, theta
    %moves(2).direction = {0.8, 1.2, pi/2}; %x, y, theta
    %moves(3).direction = {1.2, 1.6, pi}; %x, y, theta
    
    %function [minCosts, trajX, trajY, trajTht, trajV, T] = ObstacleAvoidance(trajVars, move);
%todo: add constraint for theta and omega
%plot X vs T and Y vs T separately
%label graphs correctly or use different colors
%refractor more parts: time slice duration, initial guesses

    clear all;
    MATLIBS = {'control/snopt','control/snopt/matlab/matlab', ...
               'control/optimal/Optragen/src', 'control/optimal/Optragen/', ...
               'control/optimal/'};        

    ivalab.loadLibraries(MATLIBS);
    %addpath('autogen');

    global nlp thtf gamma;
    ParamList = {'thtf', 'gamma'};

    %
    trajVars = {6, 2, 3};
    moves(1).direction = {0, 0, 0}; %x, y, theta
    moves(2).direction = {1, 1, pi/2}; %x, y, theta
    %moves(3).direction = {2.0, 1.0, pi}; %x, y, theta
    moves(3).direction = {1.5, 1.8, pi}; %x, y, theta

    hl = 1.0;
    eps = 0.0001;
    gamma = 100;

    [dim, dimSize] = size(moves);

    for i_mov = 1:dimSize-1
        %initialize variables
        x0 = moves(i_mov).direction{1};% conditions.initial{1};
        y0 = moves(i_mov).direction{2};% conditions.initial{2};
        tht0 = moves(i_mov).direction{3};% conditions.initial{3}; 

        xf = moves(i_mov+1).direction{1};% conditions.final{1};
        yf = moves(i_mov+1).direction{2};% conditions.final{2};
        thtf = moves(i_mov+1).direction{3};% conditions.final{3};

        %initial guess variables
        Time = linspace(0,1,100);
        xval = linspace(0,xf,100);
        yval = linspace(0,yf,100);
        thtval = linspace(0,thtf,100);
        vval = linspace(0,10,100);
        
        ninterv = trajVars{1};
        nsmooth = trajVars{2};
        norder = trajVars{3};

        %==[2] Define the optimization variables.
        %-- Robot state and input variables.
        %   Since theta is directly controlled, leave derivative unconstrained.
        gSym = {'x'; 'y'; 'tht';'v'};
        gpSym = {'xd'; 'yd';'thtd'}; %; 'thtd'

        %math controls and robotics

        % Create trajectory variables
        % ===========================
        x = traj(gSym{1}, ninterv,nsmooth,norder); % Arguments are ninterv, smoothness, order
        y = traj(gSym{2}, ninterv,nsmooth,norder);
        tht = traj(gSym{3}, ninterv,nsmooth,norder);
        v = traj(gSym{4}, ninterv,nsmooth,norder);

        % Create derivatives of trajectory variables
        % ==========================================
        xd  =  x.deriv(gpSym{1});
        yd  =  y.deriv(gpSym{2});
        thtd = tht.deriv(gpSym{3});

        xVars = getWorkSpaceTrajNames;

        % Define constraints
        % ==================
        Constr = constraint(x0,'x',x0,'initial', xVars) + ... % x(0)
                 constraint(y0,'y',y0,'initial', xVars) + ... % y(0)    
                 constraint(tht0,'tht',tht0,'initial', xVars) + ... %tht(0)
                 constraint(0,'v',0,'initial', xVars) + ... %v(0)
                 constraint(xf,'x',xf,'final', xVars) + ...     % x(1) = 1,  Final position, time is normalised
                 constraint(yf,'y',yf,'final', xVars) + ...     % y(1) = 1     constraint(0,'v',0,'final', xVars) + ...     % y(1) = 1
                 constraint(gamma,'gamma*cos(tht)*cos(thtf) + gamma*sin(tht)*sin(thtf)',gamma,'final', xVars) + ...                 
                 ...%constraint(eps,'tht',(2*pi),'trajectory', xVars) + ... % Dynamics as a path constraint
                 ...%constraint((-pi),'thtd',(pi),'trajectory', xVars)+ ... % Dynamics as a path constraint
                 constraint(0, 'xd - v * cos(tht)', 0, 'trajectory', xVars) + ...
                 constraint(0, 'yd - v * sin(tht)', 0, 'trajectory', xVars) + ...
                 constraint(0.1,'(x-0.4)^2 + (y-0.5)^2',Inf,'trajectory', xVars) + ... % Dynamics as a path constraint
                 constraint(0.1,'(x-1.5)^2 + (y-0.9)^2',Inf,'trajectory', xVars) + ... % Dynamics as a path constraint
                 constraint(0.1,'(x-0.8)^2 + (y-1.5)^2',Inf,'trajectory', xVars);


        % Define Cost Function
        % ====================
        Cost = cost('xd^2+yd^2','trajectory'); % Minimise energy

        % Collocation Points, using Gaussian Quadrature formula
        % =====================================================
        breaks = linspace(0,hl,ninterv+1);
        gauss = [-1 1]*sqrt(1/3)/2;
        temp = ((breaks(2:ninterv+1)+breaks(1:ninterv))/2);
        temp = temp(ones(1,length(gauss)),:) + gauss'*diff(breaks);
        colpnts = temp(:).';

        HL = [0 colpnts hl];
        %HL = linspace(0,hl,20);


        % Path where the problem related files will be stored
        % ===================================================
        pathName = './';  % Save it all in the current directory.

        % Name of the problem, will be used to identify files
        % ===================================================
        probName = 'obstacle_avoidance';

        % List of trajectories used in the problem
        % ========================================
        TrajList = traj.trajList(x,xd,y,yd,tht,thtd,v);

        nlp = ocp2nlp(TrajList, Cost,Constr, HL, ParamList,pathName,probName);
        snset('Minimize');

        xlow = -Inf*ones(nlp.nIC,1);
        xupp = Inf*ones(nlp.nIC,1);

        %Time = linspace(0,1,100);
        %xval = linspace(0,xf,100);
        %yval = linspace(0,yf,100);
        %thtval = linspace(0,thtf,100);
        %vval = linspace(0,10,100);
        xsp = createGuess(x,Time,xval);
        ysp = createGuess(y,Time,yval);
        thtsp = createGuess(tht, Time, thtval);
        vsp = createGuess(v, Time, vval);
        init = [xsp.coefs ysp.coefs thtsp.coefs vsp.coefs];% + 0.001*rand(nlp.nIC,1);
        %init = zeros(nlp.nIC,1);

        tic;
        %[x,F,inform] = snopt(init, xlow, xupp, [], [], ...
        %                     [0;nlp.LinCon.lb;nlp.nlb], [Inf;nlp.LinCon.ub;nlp.nub],...
        %                     [], [], 'ocp2nlp_cost_and_constraint');

        ghSnopt = snoptFunction(nlp);
        [x,F,inform] = snopt(init', xlow, xupp, [], [], ...
                            [0;nlp.LinCon.lb;nlp.nlb], [Inf;nlp.LinCon.ub;nlp.nub], ...
                              [], [], ...
                            ghSnopt);
                        
        toc;
        minCost = F(1);

        sp = getTrajSplines(nlp,x);
        xSP = sp{1};
        ySP = sp{2};
        thtSP = sp{3};
        vSP = sp{4};

        refinedTimeGrid = linspace(min(HL),max(HL),100);

        X = fnval(xSP,refinedTimeGrid);
        Xd = fnval(fnder(xSP),refinedTimeGrid);

        Y = fnval(ySP,refinedTimeGrid);
        Yd = fnval(fnder(ySP),refinedTimeGrid);

        Tht = fnval(thtSP, refinedTimeGrid);
        V = fnval(vSP, refinedTimeGrid);
        
        Thtd = fnval(fnder(thtSP), refinedTimeGrid);
        
        %figure(i_mov);
        %clf;
        %subplot(3,1,1); plot(X,Y);
        %text(X(100),Y(100),'\leftarrow A(\X(100),\Y(100))');
        %hold on;
        %th = 0:0.1:2*pi;
        %x = 0.4 + sqrt(.1)*cos(th);
        %y = 0.5 + sqrt(.1)*sin(th);
        %x1 = 1.5 + sqrt(.1)*cos(th);
        %y1 = 0.9 + sqrt(.1)*sin(th);
        %x2 = 0.8 + sqrt(.1)*cos(th);
        %y2 = 1.5 + sqrt(.1)*sin(th);
        %plot(x,y,'r',x1,y1,'r',x2,y2,'r');
        %axis equal;title('XY position'); xlabel('x'); ylabel('y');
        %hold off;

        %subplot(3,1,2); 
        %clf;
        %plot(Time, Tht, 'b', Time([1,end]), [thtf thtf],'r-.');hold on;
        %title('Time vs Theta'); xlabel('Time'); ylabel('Theta (rad)');
        %hold off;

        %subplot(3,1,3);
        %clf;
        %plot(Time, V);hold on;
        %title('Time vs Velocity'); xlabel('Time'); ylabel('V (m/s)');
        %hold off;

        minCosts(i_mov).val = minCost;
        trajX(i_mov).val = X;
        trajY(i_mov).val = Y;
        trajTht(i_mov).val = Tht;
        trajV(i_mov).val = V;
        T(i_mov).val = Time;
        trajThtd(i_mov).val = Thtd;
    end
    
    %plot colors or legends
    pColor = {'g'; 'b'; 'y'; 'm'; 'c'};

    %graph of X and Y trajectories
    for i_mov = 1:dimSize-1
        if i_mov == 1
            figure(i_mov);
            clf;
            th = 0:0.1:2*pi;
            x = 0.4 + sqrt(.1)*cos(th);
            y = 0.5 + sqrt(.1)*sin(th);
            x1 = 1.5 + sqrt(.1)*cos(th);
            y1 = 0.9 + sqrt(.1)*sin(th);
            x2 = 0.8 + sqrt(.1)*cos(th);
            y2 = 1.5 + sqrt(.1)*sin(th);
            plot(x,y,'r',x1,y1,'r',x2,y2,'r');
            axis equal;title('XY position'); xlabel('x'); ylabel('y');
            hold on;
        end
        X = trajX(i_mov).val; Y = trajY(i_mov).val;
        xCoord = X(end);
        yCoord = Y(end);
        plot(X,Y,pColor{i_mov});
        %text(xCoord,yCoord,'\rightarrow (\xCoord),\yCoord)');
        %hold on;       
        %hold off;
    end
    
    %graph of Time vs theta
    for i_mov = 1:dimSize-1
        if i_mov == 1
            figure(i_mov+1);
            clf;
            title('Time vs Theta'); xlabel('Time'); ylabel('Theta (rad)');
            hold on;
        end
        Time = T(i_mov).val; Tht = trajTht(i_mov).val;
        xCoord = Time(100);
        yCoord = Tht(100);
        %plot(Time, Tht, 'b', Time([1,end]), [thtf thtf],'r-.');
        plot(Time,Tht,pColor{i_mov});
        %text(xCoord,yCoord,'\leftarrow (\xCoord),\yCoord)');
        %hold on;
    end
   
    %graph of Time vs Theta
     for i_mov = 1:dimSize-1
        if i_mov == 1
            figure(i_mov+2);
            clf;
            title('Time vs Velocity'); xlabel('Time'); ylabel('V (m/s)');
            hold on;
        end
        Time = T(i_mov).val; V = trajV(i_mov).val;
        xCoord = Time(100);
        yCoord = V(100);
        plot(Time,V,pColor{i_mov});
        %text(xCoord,yCoord,'\leftarrow (\xCoord),\yCoord)');
        %hold on;
     end
    %end
    
    %graph of X vs T trajectories
    for i_mov = 1:dimSize-1
        if i_mov == 1
            figure(i_mov + 3);
            clf;
            title('Time vs X position'); xlabel('Time'); ylabel('(X)');
            hold on;
        end
        Time = T(i_mov).val; X = trajX(i_mov).val;
        xCoord = Time(end);
        yCoord = X(end);
        plot(Time,X,pColor{i_mov});
        %text(xCoord,yCoord,'\rightarrow (\xCoord),\yCoord)');
        %hold on;       
        %hold off;
    end
    
    %graph of X vs T trajectories
    for i_mov = 1:dimSize-1
        if i_mov == 1
            figure(i_mov + 4);
            clf;
            title('Time vs Y position'); xlabel('Time'); ylabel('(Y)');
            hold on;
        end
        Time = T(i_mov).val; Y = trajY(i_mov).val;
        xCoord = Time(end);
        yCoord = Y(end);
        plot(Time,Y,pColor{i_mov});
        %text(xCoord,yCoord,'\rightarrow (\xCoord),\yCoord)');
        %hold on;       
        %hold off;
    end
   
    %graph of Time vs Omega
    for i_mov=1:dimSize-1
        if i_mov == 1
            figure(i_mov +5);
            clf;
        title('Time vs Omega'); xlabel('Time'); ylabel('Omega (rad/s)');
        hold on;
        end
        Time = T(i_mov).val; Thtd =trajThtd(i_mov).val;
        plot(Time,Thtd,pColor{i_mov});        
    end
