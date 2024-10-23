function fMP = NLFD_model_parameters(MP,freq)

    if ~exist('freq','var')        
        fMP.freq = round(MP.freq(:),2);
    else
        fMP.freq = round(freq,2);
    end

    if strcmpi(MP.opt_stim,'multi')
        fMP.Fs = max(fMP.freq)*4; % sampling freq
        fMP.df = 0.1; % frequency resolution; common divisor
    elseif strcmpi(MP.opt_stim,'pure')
        fMP.Fs = fMP.freq*40;
        fMP.df = fMP.freq*0.25;
    end
    fMP.Ts = 1/fMP.Fs; % sampling period
    fMP.L = 2*ceil(fMP.Fs/fMP.df/2); % signal length 
    fMP.tt = (0:fMP.L-1)*fMP.Ts;
    fMP.dt = fMP.tt(2) - fMP.tt(1);
    fMP.ff = (0:fMP.L-1)*fMP.df; % frequency domain
    fMP.ffsym = (-fMP.L/2:fMP.L/2-1)*fMP.df; %frequency domain symmetric
    fMP.Nyqst = fMP.ff(end/2+1); %Nyquist frequency

    fMP.sind = knnsearch(fMP.ffsym.',[0;fMP.freq]).';
    fMP.freq = fMP.ffsym(fMP.sind(2:end)).';
    fMP.ind = knnsearch(fMP.ffsym.',[0;fMP.freq]).'; % setting the index to the stimulation frequency
    fMP.bins = length(fMP.ind);
end