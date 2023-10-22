% Combined Nonlinearity DF Analysis
% Authored by: Colter Russell
%
% Last updated: 08/13/2023 22:27
% 
% An analytical solution for combined amplitude and rate limit 
% nonlinearity response in the freq domain.
%
% Current Code Validity Square
%           
%               RL
%            Y     N 
%          ___________                                    
%         |     |     |                                  
%     Y   |  X  |     |                                  
% AL      |_____|_____|                                  
%         |     |     |                                  
%     N   |     |  X  |                                  
%         |_____|_____|                                  
%
% Both of the following need to be satisfied in order for the test case 
% to be valid:
%  Bif 2 occurs @ 
%    if RL/(amp*w)>1
%   Bifracation of RL sol'n occurs ~@
%       w >= (sqrt(2))*RL/amp;
%   Bifracation of AL sol'n occurs @
%       amp >= AL;


%set signal parameters 'cause we aren't amateurs!
set_case = 5;
switch set_case

case 5
    amp = 0.105;           % rad
    w   = 3;             % rad/s
case 1
    amp = 0.3;           % rad
    w   = 2;             % rad/s
case 2
    amp = 0.6;           % rad
    w   = 4;             % rad/s
case 3
    amp = 0.5;           % rad
    w   = 2;             % rad/s
case 4
    amp = 0.105;        % rad
    w   = 4;            % rad/s
otherwise
    return
end
f   = w/2/pi;       % Hz
dt  = 0.00001;      % sec
t0  = 0;            % sec
tf  = 40;
t   = t0:dt:tf;
cmd = amp*sin(w*t);
T = 1/f;

%system limits tell us what inhibits
RL  = 0.3; % Rad
AL  = 0.1; % Rad

%set loop vars (for memory we have cares)
resp      = nan(size(cmd));
colt_prev = 0;
for i = 1:length(t)

    %if limited by angle, get set to wrangle!
    colt = cmd(i); % Working var
    if abs(colt)>AL
        colt = sign(colt)*AL;
    end

    %if limited by rate, you will be late!
    if abs(colt-colt_prev)>RL*dt
        colt = colt_prev+sign(colt-colt_prev)*RL*dt;
    end

    %save the response, or rerunning will encompass your wants
    resp(i)   = colt;
    colt_prev = colt;

end

% Get a figure, plot the line
% title it up,  just in time!
% for once we have
% the figure in hand
% the grid goes on ;
% Let the calculations began!
figure
plot(t-10*T,cmd,'k-'); hold on
plot(t-10*T,resp,'g--')
title(sprintf('w = %s, Cmd = %s, AL = %s, RL = %s',num2str(w),num2str(amp),num2str(AL),num2str(RL)))
xlim([0 1/f])
grid on

% First we loop
% over each wave
% but not before 
% we floor the time array
for i = 0:floor(t(end)/f)

    % Next is part one of the magic
    % so gear up and prepare for lag, chick

    % T is the period
    % t2 is the point
    % where angle limit
    % intersects the sinisoid
    % it is used to calculate 
    % the location at t-NOT
    % we call it Y-zero
    % which brings us to the thought:
    % we need to find
    % the adjcent side length
    % also known as t1
    % we now have it all in place
    % t3 and t4
    % can be elicited via T
    % but we can also 
    % find it this fun way:

    % Order matters below...
    t2 = (i/f)+(.5/f)-asin(AL/amp)/w;
    if set_case == 4
        tstar1 = T/2 - acos(RL/(amp*w))/w;
        ystar1 = amp*sin(w*tstar1);
        tstar2 = T/2+tstar1;
        ystar2 = amp*sin(w*tstar2);
        y0 = ystar2+RL*(T-tstar2);
        t3 = tstar1+(i/f)+(ystar1+AL)/RL;
        t4 = (i/f)+T-(AL-abs(y0))/RL;
        t1 = (i/f)+(AL-y0)/RL;
    elseif set_case == 5
        tstar1 = acos(-RL/(amp*w))/w;
        %tstar1 = acos(RL/(amp*w))/w;
        ystar1 = amp*sin(w*tstar1);
        tstar2 = T/2+tstar1;
        ystar2 = amp*sin(w*tstar2);

        y0 = ystar2+RL*(T-tstar2);
        t1 = (i/f)+asin(AL/amp)/w;
        t3 = t1+T/2;
        t4 = T-asin(AL/amp)/w;

        % Testing
        %t_sq   = asin(y0/amp + sin(w*(T/2-tstar1)))/w;
        t_sq   = 2*(T-tstar2)- asin(y0);
        y_sq   = amp*sin(t_sq*w);
    else
        y0 = RL*(asin(AL/amp)/w) - AL;
        t4 = (i/f)+T-(AL-abs(y0))/RL;
        tstar1 = t2;
        ystar1 = AL;
        tstar2 = t4;
        ystar2 = -AL;
        t3 = tstar1+(i/f)+(ystar1+AL)/RL;
        t1 = (i/f)+(AL-y0)/RL;
    end


    % Add an if statement for legend prettiness
    if i ~= 1
        plot([i/f i/f],[0 y0],'bx--','HandleVisibility','off') % start of command
        plot(t1,AL,'bd','HandleVisibility','off')  %Time of switch 1
        plot(t2,AL,'bo','HandleVisibility','off')  %Time of switch 2
        plot(tstar1,ystar1,'kx','HandleVisibility','off');
        plot(t3,-AL,'b*','HandleVisibility','off') %Time of switch 3
        plot(t4,-AL,'bv','HandleVisibility','off') %Time of switch 4
        plot(tstar2,ystar2,'kx','HandleVisibility','off');
    else
        plot([i/f i/f],[0 y0],'bx--') % start of command
        plot(t1,AL,'bd')  %Time of switch 1
        plot(t2,AL,'bo')  %Time of switch 2
        plot(tstar1,ystar1,'kx');
        plot(t3,-AL,'b*') %Time of switch 3
        plot(t4,-AL,'bv') %Time of switch 4
        plot(tstar2,ystar2,'kx');
    end
    
    if set_case == 5
        plot(t_sq,y_sq,'bs','HandleVisibility','off')  %Time of switch 1
        plot(t_sq+T/2,-y_sq,'bs','HandleVisibility','off')  %Time of switch 1
    end
end

% Calculate the normal way
n_int = tf*f;
a1_out = dt*(2.*f)*trapz(resp.*cos(w.*(t)))./n_int;
b1_out = dt*(2.*f)*trapz(resp.*sin(w.*(t)))./n_int;

% Calculate Mag and Phase and Output
%mag = 20*log10(sqrt(a1_out.^2+b1_out.^2)./amp) % db
%phs = atan2(a1_out,b1_out)*180/pi % deg,
mag = sqrt(a1_out.^2+b1_out.^2)./amp % not db
phs = atan2(a1_out,b1_out) % rad

plot(t, mag*amp*sin(w*t+phs),'r--')

% If statement for not running broken code
if set_case < 4
    colt_flag = 1;
elseif set_case < 5
    colt_flag = 2;
elseif set_case == 5
    colt_flag = 3;
else
    colt_flag = 0;
end

if ~colt_flag
    legend('Command','Response','y0','t1','t2','tstar1','t3','t4','tstar2','Numerical Response')

elseif colt_flag == 3

    % Star for post intersection with AL 
    tstar1 = T/2 - acos(RL/(amp*w))/w;
    ystar1 = amp*sin(w*tstar1);
    tstar2 = T/2+tstar1;
    ystar2 = amp*sin(w*tstar2);

    % Time that AL intersects sine wave
    t2 = (T/2)-asin(AL/amp)/w;   
    % Steady state IC at t0 (o sec)
    %y0 = RL*(asin(AL/amp)/w) - AL;
    y0 = ystar2+RL*(T-tstar2);
    % Time that RL signal meets AL after t0
    t1 = (AL-y0)/RL;
    % Time that RL signal meets second AL
    %t3 = (ystar1+AL)/RL;
    t3 = tstar1+(ystar1+AL)/RL;
    % Time that AL intersects with sine wave'
    t4 = T-(AL-abs(y0))/RL;

    % Square for pre intersection with AL
    tsq1 = 2*(T-tstar2)-asin(y0);
    tsq2 = tsq1 + T/2;
    ysq1 = amp*sin(tsq1);
    ysq2 = -ysq1;

    % rewrite and pass in steps
    a1_colt = (2*RL*cos(T*w) - 2*RL + 2*RL*cos(tsq1*w) - 2*RL*cos(tsq2*w) + 2*RL*cos(tstar1*w) - 2*RL*cos(tstar2*w) + 2*w*y0*sin(tsq1*w) + 2*w*ystar1*sin(tsq2*w) - 2*w*ystar1*sin(tstar1*w) - 2*w*ystar2*sin(tstar2*w) - amp*w*cos(t1*w)^2 + amp*w*cos(t2*w)^2 - amp*w*cos(t3*w)^2 + amp*w*cos(t4*w)^2 + amp*w*cos(tsq1*w)^2 + amp*w*cos(tsq2*w)^2 - amp*w*cos(tstar1*w)^2 - amp*w*cos(tstar2*w)^2 - 2*AL*w*sin(t1*w) + 2*AL*w*sin(t2*w) + 2*AL*w*sin(t3*w) - 2*AL*w*sin(t4*w) + 2*w*ystar2*sin(T*w) + 2*RL*T*w*sin(T*w) - 2*RL*tstar2*w*sin(T*w) + 2*RL*tsq1*w*sin(tsq1*w) - 2*RL*tsq2*w*sin(tsq2*w) + 2*RL*tstar1*w*sin(tsq2*w))/(T*w^2)
    b1_colt = -(2*((RL*sin(tsq2*w) - RL*sin(tstar1*w) + w*(ystar1*cos(tsq2*w) - ystar1*cos(tstar1*w) - RL*tsq2*cos(tsq2*w) + RL*tstar1*cos(tsq2*w)))/w^2 + amp*(tsq1/2 - t1/2 + (sin(2*t1*w)/4 - sin(2*tsq1*w)/4)/w) + amp*(tsq2/2 - t3/2 + (sin(2*t3*w)/4 - sin(2*tsq2*w)/4)/w) - amp*(tstar1/2 - t2/2 + (sin(2*t2*w)/4 - sin(2*tstar1*w)/4)/w) - amp*(tstar2/2 - t4/2 + (sin(2*t4*w)/4 - sin(2*tstar2*w)/4)/w) + (w*(ystar2*cos(T*w) - ystar2*cos(tstar2*w) + RL*T*cos(T*w) - RL*tstar2*cos(T*w)) - RL*sin(T*w) + RL*sin(tstar2*w))/w^2 + (2*y0*(cos(tsq1*w)/2 - 1/2))/w - (AL*(cos(t1*w) - cos(t2*w)))/w + (AL*(cos(t3*w) - cos(t4*w)))/w - (RL*(sin(tsq1*w) - tsq1*w*cos(tsq1*w)))/w^2))/T

    % calculations
    mag_colt = sqrt(a1_colt.^2+b1_colt.^2)./amp % not db
    phs_colt = -atan2(-a1_colt,b1_colt) % rad

    % Plottem if you gottem
    plot(t, mag_colt*amp*sin(w*t+phs_colt),'r-')

    % Pull up our beautiful legend
    legend('Command','Response','y0','t1','t2','tstar1','t3','t4','tstar2','Numerical Response','Analytical Response')
elseif colt_flag == 1
    % Equations are garbage if we don't calculate:
    T = 1/f;

    % Time that AL intersects sine wave
    t2 = (T/2)-asin(AL/amp)/w;   
    % Steady state IC at t0 (o sec)
    y0 = RL*(asin(AL/amp)/w) - AL;
    % Time that RL signal meets AL after t0
    t1 = (AL-y0)/RL;
    % Time that RL signal meets second AL
    t3 = t2+2*(AL)/RL;
    % Time that AL intersects with sine wave
    t4 = T-(AL-abs(y0))/RL;

    % Half Cycle (Not as widely tested as full cycle version and Colt 
    % thinks its possible it may not include all phase degradation)
    %a1_colt  = -(4*((w*(AL*sin((T*w)/2) - AL*sin(t2*w) - (RL*T*sin((T*w)/2))/2 + RL*t2*sin((T*w)/2)) - RL*cos((T*w)/2) + RL*cos(t2*w))/w^2 - (RL*(2*sin((t1*w)/2)^2 - t1*w*sin(t1*w)))/w^2 + (y0*sin(t1*w))/w - (AL*(sin(t1*w) - sin(t2*w)))/w))/T
    %b1_colt  = -(4*((RL*sin((T*w)/2) - RL*sin(t2*w) + w*(AL*cos((T*w)/2) - AL*cos(t2*w) - (RL*T*cos((T*w)/2))/2 + RL*t2*cos((T*w)/2)))/w^2 + (2*y0*(cos(t1*w)/2 - 1/2))/w - (AL*(cos(t1*w) - cos(t2*w)))/w - (RL*(sin(t1*w) - t1*w*cos(t1*w)))/w^2))/T

    % Full cycle version (likely can be reduced)
    %a1_colt = -(2*((RL*cos(t3*w) - RL*cos(t2*w) + w*(AL*sin(t2*w) - AL*sin(t3*w) - RL*t2*sin(t3*w) + RL*t3*sin(t3*w)))/w^2 + (w*(AL*sin(T*w) - AL*sin(t4*w) - RL*T*sin(T*w) + RL*t4*sin(T*w)) - RL*cos(T*w) + RL*cos(t4*w))/w^2 + (RL*(2*sin((t1*w)/2)^2 - t1*w*sin(t1*w)))/w^2 - (y0*sin(t1*w))/w + (AL*(sin(t1*w) - sin(t2*w)))/w - (AL*(sin(t3*w) - sin(t4*w)))/w))/T;
    %b1_colt = (2*((RL*sin(T*w) - RL*sin(t4*w) + w*(AL*cos(T*w) - AL*cos(t4*w) - RL*T*cos(T*w) + RL*t4*cos(T*w)))/w^2 + (RL*sin(t2*w) - RL*sin(t3*w) + w*(AL*cos(t2*w) - AL*cos(t3*w) - RL*t2*cos(t3*w) + RL*t3*cos(t3*w)))/w^2 - (2*y0*(cos(t1*w)/2 - 1/2))/w + (AL*(cos(t1*w) - cos(t2*w)))/w - (AL*(cos(t3*w) - cos(t4*w)))/w + (RL*(sin(t1*w) - t1*w*cos(t1*w)))/w^2))/T;

    % Simplified expression
    a1_colt = (4*RL/T/w^2)*(cos(t1*w) + cos(t2*w));
    b1_colt = (4*RL/T/w^2)*(sin(t1*w) + sin(t2*w));

    % calculations
    mag_colt = sqrt(a1_colt.^2+b1_colt.^2)./amp % not db
    phs_colt = -atan2(-a1_colt,b1_colt) % rad

    % Plottem if you gottem
    plot(t, mag_colt*amp*sin(w*t+phs_colt),'r-')

    % Pull up our beautiful legend
    legend('Command','Response','y0','t1','t2','tstar1','t3','t4','tstar2','Numerical Response','Analytical Response')
elseif colt_flag == 2

    tstar1 = T/2 - acos(RL/(amp*w))/w;
    ystar1 = amp*sin(w*tstar1);
    tstar2 = T/2+tstar1;
    ystar2 = amp*sin(w*tstar2);

    % Time that AL intersects sine wave
    t2 = (T/2)-asin(AL/amp)/w;   
    % Steady state IC at t0 (o sec)
    %y0 = RL*(asin(AL/amp)/w) - AL;
    y0 = ystar2+RL*(T-tstar2);
    % Time that RL signal meets AL after t0
    t1 = (AL-y0)/RL;
    % Time that RL signal meets second AL
    %t3 = (ystar1+AL)/RL;
    t3 = tstar1+(ystar1+AL)/RL;
    % Time that AL intersects with sine wave'
    t4 = T-(AL-abs(y0))/RL;

    % Full cycle version w/ tstar
    %a1_colt = (2*((cos(2*t2*w) - cos(2*tstar1*w))/(4*w) + (cos(2*t4*w) - cos(2*tstar2*w))/(4*w) + (w*(ystar1*sin(t3*w) - ystar1*sin(tstar1*w) - RL*t3*sin(t3*w) + RL*tstar1*sin(t3*w)) - RL*cos(t3*w) + RL*cos(tstar1*w))/w^2 + (RL*cos(T*w) - RL*cos(tstar2*w) + w*(ystar2*sin(T*w) - ystar2*sin(tstar2*w) + RL*T*sin(T*w) - RL*tstar2*sin(T*w)))/w^2 - (RL*(2*sin((t1*w)/2)^2 - t1*w*sin(t1*w)))/w^2 + (y0*sin(t1*w))/w - (AL*(sin(t1*w) - sin(t2*w)))/w + (AL*(sin(t3*w) - sin(t4*w)))/w))/T
    %b1_colt = -(2*(t2/2 + t4/2 - tstar1/2 - tstar2/2 + (RL*sin(t3*w) - RL*sin(tstar1*w) + w*(ystar1*cos(t3*w) - ystar1*cos(tstar1*w) - RL*t3*cos(t3*w) + RL*tstar1*cos(t3*w)))/w^2 + (w*(ystar2*cos(T*w) - ystar2*cos(tstar2*w) + RL*T*cos(T*w) - RL*tstar2*cos(T*w)) - RL*sin(T*w) + RL*sin(tstar2*w))/w^2 - (sin(2*t2*w)/4 - sin(2*tstar1*w)/4)/w - (sin(2*t4*w)/4 - sin(2*tstar2*w)/4)/w + (2*y0*(cos(t1*w)/2 - 1/2))/w - (AL*(cos(t1*w) - cos(t2*w)))/w + (AL*(cos(t3*w) - cos(t4*w)))/w - (RL*(sin(t1*w) - t1*w*cos(t1*w)))/w^2))/T

% rewrite and pass in steps
    a1_colt = (2*((w*(ystar1*sin(t3*w) - ystar1*sin(tstar1*w) - RL*t3*sin(t3*w) + RL*tstar1*sin(t3*w)) - RL*cos(t3*w) + RL*cos(tstar1*w))/w^2 + (RL*cos(T*w) - RL*cos(tstar2*w) + w*(ystar2*sin(T*w) - ystar2*sin(tstar2*w) + RL*T*sin(T*w) - RL*tstar2*sin(T*w)))/w^2 - (RL*(2*sin((t1*w)/2)^2 - t1*w*sin(t1*w)))/w^2 + (y0*sin(t1*w))/w - (AL*(sin(t1*w) - sin(t2*w)))/w + (AL*(sin(t3*w) - sin(t4*w)))/w + (amp*(cos(2*t2*w) - cos(2*tstar1*w)))/(4*w) + (amp*(cos(2*t4*w) - cos(2*tstar2*w)))/(4*w)))/T
    b1_colt = -(2*((RL*sin(t3*w) - RL*sin(tstar1*w) + w*(ystar1*cos(t3*w) - ystar1*cos(tstar1*w) - RL*t3*cos(t3*w) + RL*tstar1*cos(t3*w)))/w^2 - amp*(tstar1/2 - t2/2 + (sin(2*t2*w)/4 - sin(2*tstar1*w)/4)/w) - amp*(tstar2/2 - t4/2 + (sin(2*t4*w)/4 - sin(2*tstar2*w)/4)/w) + (w*(ystar2*cos(T*w) - ystar2*cos(tstar2*w) + RL*T*cos(T*w) - RL*tstar2*cos(T*w)) - RL*sin(T*w) + RL*sin(tstar2*w))/w^2 + (2*y0*(cos(t1*w)/2 - 1/2))/w - (AL*(cos(t1*w) - cos(t2*w)))/w + (AL*(cos(t3*w) - cos(t4*w)))/w - (RL*(sin(t1*w) - t1*w*cos(t1*w)))/w^2))/T

    % calculations
    mag_colt = sqrt(a1_colt.^2+b1_colt.^2)./amp % not db
    phs_colt = -atan2(-a1_colt,b1_colt) % rad

    % Plottem if you gottem
    plot(t, mag_colt*amp*sin(w*t+phs_colt),'r-')

    % Pull up our beautiful legend
    legend('Command','Response','y0','t1','t2','tstar1','t3','t4','tstar2','Numerical Response','Analytical Response')
end
