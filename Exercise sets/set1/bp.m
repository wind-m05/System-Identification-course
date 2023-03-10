function [mag,phase_var,freq] = bp(data,options,range)
%BP is a custom bodeplotter
datatypetest = class(data);
if strcmp(datatypetest,'tf')
data = frd(data,range);
end
vec = data.responsedata;
mag = squeeze(abs(vec));
phase_var = rad2deg(squeeze(unwrap(angle(vec))));
freq = data.frequency; 
% if options.unwrap == true
% phase_var = unwrap(phase_var);
% end

% Plotting
if options.subplot == true % Then do mag and phase
    % Scaling of plot
    if strcmp(options.scale,'mag')
        if strcmp(options.plot, 'scatter')
            subplot(2,1,1)
            scatter(freq,mag)
            set(gca,'xscale','log')
            subplot(2,1,2)
            scatter(freq,phase_var)
            set(gca,'xscale','log')
        else % Else make a standard bode plot
            subplot(2,1,1)
            semilogx(freq,mag)
            subplot(2,1,2)
            semilogx(freq,phase_var)
        end
    else % else do log and db
        if strcmp(options.plot, 'scatter')
            subplot(2,1,1)
            scatter(freq,mag2db(mag))
            set(gca,'xscale','log')
            subplot(2,1,2)
            scatter(freq,phase_var)
            set(gca,'xscale','log')
        else % Else make a standard bode plot
            subplot(2,1,1)
            semilogx(freq,mag2db(mag))
            subplot(2,1,2)
            semilogx(freq,phase_var)
        end 
    end    
    
else 
    % Scaling of plot
    if strcmp(options.scale,'mag')
        if strcmp(options.plot, 'scatter')
            scatter(freq,mag)
            set(gca,'xscale','log')
        else % Else make a standard bode plot
            semilogx(freq,mag)
        end
    else % else do log and db
        if strcmp(options.plot, 'scatter')
            scatter(freq,mag2db(mag))
            set(gca,'xscale','log')
        else % Else make a standard bode plot
            semilogx(freq,mag2db(mag))
        end 
    end    
end
end

