function [Nd,El,MP,fMP,FLD,OHC,R] = NLFD_virtualCochleaHarmonic(z)

%     set_Global_coefficients; % this is setting coefficients for parameter study

    To = onset_msg;

    opt_save = 0;

    fftw('dwisdom',[]);
    fftw('planner','patient');

    pool = gcp('nocreate');
    workers = 4;
    if ~isempty(pool)
        if ne(pool.NumWorkers, workers)
            delete(pool);
            parpool('local',workers);
        end
    else
        parpool('local',workers);
    end   

    global coef  

    if ~exist('z','var')
        z = '';
    end
    MP = set_MP(z);
    MP.coef = coef; % save coef information before save

    [Nd,El,OHC,IHC,MP,FLD] = initiate(MP);

    PoI = define_PoI(MP,Nd,El);

    % frequency independent part. only contains inviscid scala fluid formulation
    [A0,AMat,El,FLD,OHC,Ce,Ge,Ie,V0] = assemble_A0(Nd,El,OHC,IHC,FLD,MP);
    
    if strcmp(MP.opt_stim,'pure') 
        [Nd,El,MP,fMP,FLD,OHC,R] = NLFD_pure_tone(MP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0);
    elseif strcmp(MP.opt_stim,'multi') 
        [Nd,El,MP,fMP,FLD,OHC,R] = NLFD_multi_tone(MP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0);
    end

    if opt_save
        ct = clock;
        date_str = [num2str(ct(2),'%02.2d'),num2str(ct(3),'%02.2d'),num2str(rem(ct(1),2000),'%02.2d')];
        time_str = [num2str(ct(4),'%02.2d'),num2str(ct(5),'%02.2d')];
        formatSpec = '%.2f';
        file_name = ['./VC_code/houtput/',date_str,'_',time_str,'.mat'];
        save(file_name,'MP','fMP','Nd','El','OHC','FLD','R');
    end

    fin_msg(To);

end

function [Nd,El,MP,fMP,FLD,OHC,R] = NLFD_pure_tone(MP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0)

    ct = clock;

    R.freq = reshape(MP.freq,1,[]);   % initiate R

    for i = 1:MP.Nf

        fMP = NLFD_model_parameters(MP,MP.freq(i));

        Uf = zeros(MP.tdof,fMP.L,'like',1j);        
        Uf(end) = 1j*1e-13;
        ind = (fMP.ffsym == 0);
        idx = MP.dof(5).odof;
        HB = [OHC.HB];
        Uf(idx,ind) = [HB.po];
        idx = MP.dof(3).edof;
        Uf(idx,ind) = V0;

        MP.opt_structure_relaxation = 0;
        MP.opt_fluid_relaxation = 0;
        MP.key_fOHC = 0; MP.key_fMET = 0;
        MP.Estim = 0;
        opt_lin = 1; % 1: linear; 0: nonlin 
        [Uf,AMat] = NLFD_model(MP,fMP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0,Uf,opt_lin);
        R_i = store_results('NLFD_psv',ct,PoI,MP,fMP,Nd,El,OHC,FLD,AMat,Uf);

        idx = find(fMP.ffsym(fMP.ind) == fMP.freq);
        fields = fieldnames(R_i);
        for f = 1:length(fields)
            R(1).stim = 'psv';
            R(1).(fields{f})(:,i) = R_i.(fields{f})(:,idx);
        end

        if MP.sensitive_cochlea
            MP.key_fOHC = 1; MP.key_fMET = 1;
            if MP.vMC == 1, MP.Estim = 1; end
            [Uf,AMat] = NLFD_model(MP,fMP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0,Uf,1);
            if MP.nonlin
                MP.opt_structure_relaxation = 1;
                MP.opt_fluid_relaxation = 1;
                [Uf,AMat] = NLFD_model(MP,fMP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0,Uf,0);
            end                        
            R_i = store_results('NLFD_act',ct,PoI,MP,fMP,Nd,El,OHC,FLD,AMat,Uf);
            idx = find(fMP.ffsym(fMP.ind) == fMP.freq);
            fields = fieldnames(R_i);
            for f = 1:length(fields)
                R(2).stim = 'act';
                R(2).(fields{f})(:,i) = R_i.(fields{f})(:,idx);
            end
        end

    end

end

function [Nd,El,MP,fMP,FLD,OHC,R] = NLFD_multi_tone(MP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0)

    ct = clock;

    fMP = NLFD_model_parameters(MP);

    Uf = zeros(MP.tdof,fMP.L,'like',1j);        
    Uf(end) = 1j*1e-13;
    ind = (fMP.ffsym == 0);
    idx = MP.dof(5).odof;
    HB = [OHC.HB];
    Uf(idx,ind) = [HB.po];
    idx = MP.dof(3).edof;
    Uf(idx,ind) = V0;

    MP.opt_structure_relaxation = 0;
    MP.opt_fluid_relaxation = 0;
    MP.key_fOHC = 0; MP.key_fMET = 0;
    [Uf,AMat] = NLFD_model(MP,fMP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0,Uf);
    R = store_results('NLFD_psv',ct,PoI,MP,fMP,Nd,El,OHC,FLD,AMat,Uf,0);

    if MP.sensitive_cochlea
        MP.opt_structure_relaxation = 0;
        MP.opt_fluid_relaxation = 0;
        MP.key_fOHC = 1; MP.key_fMET = 1;
        [Uf,AMat] = NLFD_model(MP,fMP,PoI,Nd,El,OHC,FLD,A0,AMat,Ce,Ge,Ie,V0,Uf);
        store_results('NLFD_act',ct,PoI,MP,fMP,Nd,El,OHC,FLD,AMat,Uf,1);
    end

end

function To = onset_msg
To=clock;
fprintf(1,'\n\n   Performing organ of Corti FE Harmonic Analysis...\n\n');
fprintf(1,'  Program starts at %d:%d\n\n',To(4:5));

end % of function onset_msg()

function fin_msg(To)

Tf = clock;
eT=etime(Tf,To);
secs=rem(round(eT),60);
mins=rem(floor(eT/60),60);
hrs=floor(eT/3600);

fprintf(1,'  Program ends at %3d:%3d:%3d sec \n',round(Tf(4:6)));
fprintf(1,'  It took %d hr %d min %d sec total.',hrs,mins,secs);
fprintf(1,'\n=======================================================\n');

end % of function fin_msg()
