function [] = plot_mvdr( name )
load(name)
    numst = 1000; % resolution in spatial response
%    W_H = conj(W(Ndata, :)); % Hemitian transpose of last one
    st = linspace(-1,1,numst); % sine(theta) space
    est = exp(-1j*pi*[0:(5-1)]'*st); % steering matrix
    % amplitude response
    P = 20*log10(abs(W_H*est).^2); 
    plot(st,P)
end

