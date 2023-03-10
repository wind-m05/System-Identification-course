function [mag,phase_var,freq,l1,l2,p1,p2] = bp(data,options,range)
%BP is a custom bodeplotter
datatypetest = class(data);
if strcmp(datatypetest,'tf')
data = frd(data,range);
end
vec = data.responsedata;
mag = squeeze(abs(vec));
phase_var = rad2deg(squeeze(unwrap(angle(vec))));
freq = data.frequency; 

% Plotting
if options.subplot == true % Then do mag and phase
    % Scaling of plot
    if strcmp(options.yscale,'mag')
        if strcmp(options.plot, 'scatter')
            p1 = subplot(2,1,1);
            l1 = scatter(freq,mag);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
            p2 = subplot(2,1,2);
            l2 = scatter(freq,phase_var);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        else % Else make a standard bode plot
            p1 = subplot(2,1,1);
            l1 = plot(freq,mag);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
            p2 = subplot(2,1,2);
            l2 = plot(freq,phase_var);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        end
    else % else do log and db
        if strcmp(options.plot, 'scatter')
            p1 = subplot(2,1,1);
            l1 = scatter(freq,mag2db(mag));
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
            p2 = subplot(2,1,2);
            l2 = scatter(freq,phase_var);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        else % Else make a standard bode plot
            p1 = subplot(2,1,1);
            l1 = plot(freq,mag2db(mag));
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
            p2 = subplot(2,1,2);
            l2 = plot(freq,phase_var);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        end 
    end    
else 
    % Scaling of plot
    if strcmp(options.yscale,'mag')
        if strcmp(options.plot, 'scatter')
            l1 = scatter(freq,mag);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        else % Else make a standard bode plot
            l1 = plot(freq,mag);
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        end
    else % else do log and db
        if strcmp(options.plot, 'scatter')
            l1 = scatter(freq,mag2db(mag));
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        else % Else make a standard bode plot
            l1 = plot(freq,mag2db(mag));
            if strcmp(options.xscale, 'log')
            set(gca,'xscale','log')
            end
        end 
    end    
end
end

