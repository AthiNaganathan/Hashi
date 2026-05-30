function FesCalc_Block_gen(pname, struct_file, ene, DS, DCp, T, IS)
    %% Looping PDB files
    
    for run=1    
        pdb=char(pname); % input PDB
        stridefile = char(struct_file); % input struct.txt from STRIDE
        % DS=-14.5/1000; % in kJ/mol/K... entropic penalty
        % ene=-79/1000; % kJ/mol... van der Waals interaction energy per native contact
        % DCp=-0.36/1000; % kJ/mol/K per native contact... heat capacity change per native contact
        % disr = [1,4:5,10:16,20:23,25:26,32:33,38:44,58:59,64:78,84,90:91,96:106,111:113,120:121,125:128,...
        %         133:135,145:148,152:158,160:162,177,179:182,190:191,198,199,203:206,208,213:215,226:235,237:240,242:245,247:250,252:255,257:260,263,268:269,...
        %         275:276,280,282:284,290:293,299:300,307,313:319,321:323,331:333,337,344,345,347,354:357,361,365:368,370:371,373,...
        %         376:377,385:387,391:395,403:404,409:411,415,421:423,428:429,431:437,446:448,459,462,...
        %         471:472,481:484,486:490,496,498:502,507:508,511,516:517,525:528,533:534,537]; % residues that are either coil or glycine
        % ppos = [24,50,159,163,164,178,183,236,241,246,251,256,281,298,305,320,324,346,348,369,372,378,430,485,497,538,539]; % positions of proline residues
        % T=298; % Temperature in K 
        % IS=0.05; % Ionic strength in Molar units
        C=0; % denaturant concentration, not used and kept fixed at zero.
        disr = []; % residues that are either coil or glycine
        ppos = []; % positions of proline residues  
    
        %% Contact Maps obtained through PDB

        aa2={'GLY','ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP','SER','ASP','ASN','THR','GLU','GLN','HIS','LYS','ARG','PRO','CYS'};
        aacode={'G','A','V','L','I','M','F','Y','W','S','D','N','T','E','Q','H','K','R','P','C'};
        stridefile=fullfile(pwd, '..', 'data/', stridefile);
        if ~isfile(stridefile)
            error("ERROR : File %s not found.", stridefile)
        end
        hh9=fopen(stridefile(1,:),'rt');
        lkk=fgetl(hh9);
        lengthprot = 1;
        while lkk > 0
            if length(lkk)> 25 && strcmp(lkk(1:3),'ASG')
                if strcmp(lkk(6:8),'PRO')
                    ppos(end+1)=str2double(lkk(17:20));
                elseif strcmp(lkk(6:8),'GLY') || (strcmp(lkk(25),'H')==0 && strcmp(lkk(25),'E')==0 && strcmp(lkk(25),'G')==0)
                    disr(end+1)=str2double(lkk(17:20));
                end
                [~,x2]=find(strcmp(lkk(6:8),aa2));
                protseq(1,lengthprot)=aacode(x2);
                lengthprot = lengthprot+1;
            end
            lkk=fgetl(hh9);
        end
        fclose(hh9);

        aa=pdb(1,:);

        opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', ['contactmapmatElecB',aa,'.dat']);
        load(opPath)
        eval(['cmapmask=contactmapmatElecB',aa,';']);

        opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', ['contdistElecB',aa,'.dat']);
        load(opPath)
        eval(['contdist=contdistElecB',aa,';']);
        
        opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', [aa '_BlockDet.dat']);
        load(opPath)
        eval(['BlockDet=',aa,'_BlockDet;']);

        opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', ['BlockSize',aa,'.dat']);
        load(opPath)
        eval(['bs=BlockSize',aa,';']);

        nres=length(cmapmask);
        nres1 = length(unique(BlockDet(:,1)));
        
        %% Variables
        
        wi=4;             % for number of calculating windows
        Mw=nres1*110;
        
        
        %% Constants
        R=0.008314;
        Tref=385;
        zval=exp(DS./R);
        zvalc=exp((DS-(6.0606/1000))./R);
        zjj=zval.*ones(nres1,1);
        zjj(disr)=zvalc;
        zjj(ppos)=1;
        zvec = ones(BlockDet(end,2),1);
        for i=1:length(BlockDet)
            zvec(BlockDet(i,2)) = zvec(BlockDet(i,2)) * zjj(BlockDet(i,1));
        end
        
        %% calculate Zfin
        Zfin=zeros(length(T),1);
        for kk=1:length(T)
            for ll=1:length(C)
                %% Elec. energy involving IS contribution
                ISfac=5.66*sqrt(IS/T(kk))*sqrt(80/29);
                emapmask=zeros(nres);
                for i=1:nres
                    for iin=i:nres
                        x1=find(contdist(:,1)==i & contdist(:,2)==iin);
                        emapmask(i,iin)=sum(contdist(x1,5).*exp(-ISfac.*contdist(x1,3)));
                    end
                end
                pepval=zeros(nchoosek(nres+1,2)+2*nchoosek(nres+1,4),7);
                k=1;
                %% Generating Combinations for SSA
                for i=1:nres
                    for iin=1:nres-i+1
                        sw=0;
                        stabEtemp=0; eneEtemp=0; nconttemp=0;
                        wii=floor(iin/wi);
                        swgen=zeros(wi,1);
                        wistart=i;
                        for wiii=1:wi-1
                            nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,i:i+iin-1)));
                            stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                            eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,i:i+iin-1)));
                            swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                            wistart=i+wii*wiii+1;
                        end
                        nconttemp=sum(sum(cmapmask(wistart:i+iin-1,i:i+iin-1)));
                        stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                        eneEtemp=sum(sum(emapmask(wistart:i+iin-1,i:i+iin-1)));
                        swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                        sw=prod(swgen);
                        
                        pepval(k,:)=[(iin) sw i iin 0 0 1];
                        k=k+1;
                    end
                    fprintf("%d ; %d ; %d ; %d (a) \n", run, kk, ll, i);
                end
                
                %% Generating Combinations for DSA
                for i=1:nres
                    for iin=1:nres-i+1
                        for j=i+iin+1:nres
                            for jin=1:nres-j+1
                                
                                %for island 1
                                sw1=0;
                                stabEtemp=0; eneEtemp=0; nconttemp=0;
                                wii=floor(iin/wi);
                                swgen=zeros(wi,1);
                                wistart=i;
                                for wiii=1:wi-1
                                    nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,i:i+iin-1)));
                                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                    eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,i:i+iin-1)));
                                    swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                                    wistart=i+wii*wiii+1;
                                    stabEtemp=0; eneEtemp=0; nconttemp=0;
                                end
                                nconttemp=sum(sum(cmapmask(wistart:i+iin-1,i:i+iin-1)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:i+iin-1,i:i+iin-1)));
                                swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                                sw1=prod(swgen);
                                
                                %for island 2
                                sw2=0;
                                stabEtemp=0; eneEtemp=0; nconttemp=0;
                                wii=floor(jin/wi);
                                swgen=zeros(wi,1);
                                wistart=j;
                                for wiii=1:wi-1
                                    nconttemp=sum(sum(cmapmask(wistart:j+wii*wiii,j:j+jin-1)));
                                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                    eneEtemp=sum(sum(emapmask(wistart:j+wii*wiii,j:j+jin-1)));
                                    swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+wii*wiii));
                                    wistart=j+wii*wiii+1;
                                end
                                nconttemp=sum(sum(cmapmask(wistart:j+jin-1,j:j+jin-1)));
                                stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneEtemp=sum(sum(emapmask(wistart:j+jin-1,j:j+jin-1)));
                                swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+jin-1));
                                sw2=prod(swgen);
                                
                                sw=sw1.*sw2;
                                
                                pepval(k,:)=[(iin+jin) sw i iin j jin 2];
                                k=k+1;
                            end
                        end
                    end
                    fprintf("%d ; %d ; %d ; %d (b) \n", run, kk, ll, i);
                end
                
                %% Generating Combinations for DSAw/L
                for i=1:nres
                    for iin=1:nres-i+1
                        for j=i+iin+1:nres
                            for jin=1:nres-j+1
                                stabE=0; eneE=0; sw=0;
                                vv = [(i:i+iin-1) (j:j+jin-1)];
                                ncont=sum(sum(cmapmask(vv,vv)));
                                stabE=(ncont*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                eneE=sum(sum(emapmask(vv,vv)));
                                if sum(sum(cmapmask(i:i+iin-1,j:j+jin-1)))~=0 || sum(sum(emapmask(i:i+iin-1,j:j+jin-1)))~=0
                                    %from island 1
                                    sw1=0;
                                    stabEtemp=0; eneEtemp=0; nconttemp=0;
                                    wii=floor(iin/wi);
                                    swgen=zeros(wi,1);
                                    wistart=i;
                                    for wiii=1:wi-1
                                        nconttemp=sum(sum(cmapmask(wistart:i+wii*wiii,vv)));
                                        stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                        eneEtemp=sum(sum(emapmask(wistart:i+wii*wiii,vv)));
                                        swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+wii*wiii));
                                        wistart=i+wii*wiii+1;
                                    end
                                    nconttemp=sum(sum(cmapmask(wistart:i+iin-1,vv)));
                                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                    eneEtemp=sum(sum(emapmask(wistart:i+iin-1,vv)));
                                    swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:i+iin-1));
                                    sw1=prod(swgen);
                                    
                                    %from island 2
                                    sw2=0;
                                    stabEtemp=0; eneEtemp=0; nconttemp=0;
                                    wii=floor(jin/wi);
                                    swgen=zeros(wi,1);
                                    wistart=j;
                                    for wiii=1:wi-1
                                        nconttemp=sum(sum(cmapmask(wistart:j+wii*wiii,vv)));
                                        stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                        eneEtemp=sum(sum(emapmask(wistart:j+wii*wiii,vv)));
                                        swgen(wiii,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+wii*wiii));
                                        wistart=j+wii*wiii+1;
                                    end
                                    nconttemp=sum(sum(cmapmask(wistart:j+jin-1,vv)));
                                    stabEtemp=(nconttemp*(ene+DCp*(T(kk)-Tref)-T(kk)*DCp*log(T(kk)/Tref)));
                                    eneEtemp=sum(sum(emapmask(wistart:j+jin-1,vv)));
                                    swgen(wi,1)=exp(-(stabEtemp+eneEtemp)/(R*T(kk)))*prod(zvec(wistart:j+jin-1));
                                    sw2=prod(swgen);
                                    
                                    sw=sw1*sw2*zvalc^(j-(i+iin));
                                else
                                    sw = 0;
                                end
                                pepval(k,:)=[(iin+jin) sw i iin j jin 3];
                                k=k+1;
                            end
                        end
                    end
                    fprintf("%d ; %d ; %d ; %d (c) \n", run, kk, ll, i);
                end
                
                Zfin(kk) = sum(pepval(:,2));
            end
        end

        pepval = pepval(pepval(:,2)~=0,:);
        %% 1D free energy profile
        fes = zeros(nres,1);
        for i=1:nres
            fes(i) = sum(pepval(pepval(:,1)==i,2));
        end
        fes = fes./sum(fes);
        fes = -R*T*log(fes);
        
        %% Major Folding Path (Residue folding probability as the function of RC) and 2D free-energy profile
        hv = round(nres/2);
        Fpath = zeros(nres1,nres);
        FpathZ = zeros(nres,1);
        conv = zeros(nres,2);
        fes2D = zeros(hv+1,nres-hv+1);
        for i=1:nres
            f = find(BlockDet(:,2)==i);
            conv(i,:) = [BlockDet(f(1),1) BlockDet(f(end),1)];
        end
        
        for i=1:length(pepval)
            pept = zeros(nres1,1);
            pept(conv(pepval(i,3),1):conv(pepval(i,3)+pepval(i,4)-1,2))=1;
            if pepval(i,end)>1; pept(conv(pepval(i,5),1):conv(pepval(i,5)+pepval(i,6)-1,2))=1; end
            Fpath(:,pepval(i,1)) = Fpath(:,pepval(i,1)) + pepval(i,2)*pept;
            FpathZ(pepval(i,1)) = FpathZ(pepval(i,1)) + pepval(i,2);
            pept = zeros(nres,1);
            pept(pepval(i,3):pepval(i,3)+pepval(i,4)-1)=1;
            if pepval(i,end)>1; pept(pepval(i,5):pepval(i,5)+pepval(i,6)-1)=1; end
            fes2D(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1) = fes2D(sum(pept(1:hv))+1,sum(pept(hv+1:end))+1)+pepval(i,2);
        end
        for i=1:nres
            Fpath(:,i) = Fpath(:,i)./FpathZ(i);
        end
        
        %% 2D free energy profile
        fes2D = fes2D./sum(sum(fes2D));
        fes2D = -R*T*log(fes2D);
    end

    fes = [(1:nres)' fes];
    
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', '1D_FreeEnergyProfile.txt');
    ff = fopen(opPath,'wt');
    fprintf(ff,'REMARK @ Input Parameters\n');
    fprintf(ff,'REMARK @     van der Waals interaction energy (J/mol per native contact): %3.2f\n',ene*1000);
    fprintf(ff,'REMARK @     Entropic Cost (J/mol.K per residue)                        : %3.2f\n',DS*1000);
    fprintf(ff,'REMARK @     Heat Capacity Change (J/mol.K per native contact)          : %3.2f\n',DCp*1000);
    fprintf(ff,'REMARK @     Temperature (K)                                            : %d\n',T);
    fprintf(ff,'REMARK @     Ionic Strenght (M)                                         : %1.2f\n',IS);
    fprintf(ff,'REMARK @     Block Size                                                 : %d\n',bs);
    fprintf(ff,'REMARK %% Title "1D Free Energy Profile"\n');
    fprintf(ff,'REMARK %% Legend \n');
    fprintf(ff,'REMARK %%    column 1 "No. of Structured Blocks"  \n');
    fprintf(ff,'REMARK %%    column 2 "Free Energy (kJ/mol)"  \n');
    for i=1:length(fes)
        fprintf(ff,'%3d %3.2f\n',fes(i,1),fes(i,2));
    end
    fclose(ff);
    
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', '2D_FreeEnergySurface.txt');
    ff = fopen(opPath,'wt');
    fprintf(ff,'REMARK @ Input Parameters\n');
    fprintf(ff,'REMARK @     van der Waals interaction energy (J/mol per native contact): %3.2f\n',ene*1000);
    fprintf(ff,'REMARK @     Entropic Cost (J/mol.K per residue)                        : %3.2f\n',DS*1000);
    fprintf(ff,'REMARK @     Heat Capacity Change (J/mol.K per native contact)          : %3.2f\n',DCp*1000);
    fprintf(ff,'REMARK @     Temperature (K)                                            : %d\n',T);
    fprintf(ff,'REMARK @     Ionic Strenght (M)                                         : %1.2f\n',IS);
    fprintf(ff,'REMARK @     Block Size                                                 : %d\n',bs);
    fprintf(ff,'REMARK %% Title "2D Free Energy Surface"\n');
    fprintf(ff,'REMARK %% Legend \n');
    fprintf(ff,'REMARK %%    column 1 "No. of Structured Blocks in the N-termini"  \n');
    fprintf(ff,'REMARK %%    column 2 "No. of Structured Blocks in the C-termini"  \n');
    fprintf(ff,'REMARK %%    column 3 "Free Energy (kJ/mol)"  \n');
    fprintf(ff,'REMARK %%    Use Plot_2D_FE_Surface.m to plot figure as shown in the pFolMech web server  \n');
    for i=1:length(fes2D(:,1))
        for j=1:length(fes2D(1,:))
            fprintf(ff,'%3d %3d %3.2f\n',i-1,j-1,fes2D(i,j));
        end
    end
    fclose(ff);
    
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'ResFoldProb_vs_RC.txt');
    ff = fopen(opPath,'wt');
    fprintf(ff,'REMARK @ Input Parameters\n');
    fprintf(ff,'REMARK @     van der Waals interaction energy (J/mol per native contact): %f\n',ene*1000);
    fprintf(ff,'REMARK @     Entropic Cost (J/mol.K per residue)                        : %3.2f\n',DS*1000);
    fprintf(ff,'REMARK @     Heat Capacity Change (J/mol.K per native contact)          : %3.2f\n',DCp*1000);
    fprintf(ff,'REMARK @     Temperature (K)                                            : %d\n',T);
    fprintf(ff,'REMARK @     Ionic Strenght (M)                                         : %1.2f\n',IS);
    fprintf(ff,'REMARK @     Block Size                                                 : %d\n',bs);
    fprintf(ff,'REMARK %% Title "Folding Probability as the function of RC"\n');
    fprintf(ff,'REMARK %% Legend \n');
    fprintf(ff,'REMARK %%    column 1 "No. of Structured Blocks"  \n');
    fprintf(ff,'REMARK %%    column 2 "Residue Index"  \n');
    fprintf(ff,'REMARK %%    column 3 "Folding Probability"  \n');
    fprintf(ff,'REMARK %%    Use Plot_RCvsResFoldProb.m to plot figure as shown in the pFolMech web server  \n');
    for i=1:length(Fpath(:,1))
        for j=1:length(Fpath(1,:))
            fprintf(ff,'%3d %3d %1.3f\n',j-1,i,Fpath(i,j));
        end
    end
    fclose(ff);
    
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'Microstates.txt');
    ff = fopen(opPath,'wt');
    fprintf(ff,'REMARK @ Input Parameters\n');
    fprintf(ff,'REMARK @     van der Waals interaction energy (J/mol per native contact): %f\n',ene*1000);
    fprintf(ff,'REMARK @     Entropic Cost (J/mol.K per residue)                        : %3.2f\n',DS*1000);
    fprintf(ff,'REMARK @     Heat Capacity Change (J/mol.K per native contact)          : %3.2f\n',DCp*1000);
    fprintf(ff,'REMARK @     Temperature (K)                                            : %d\n',T);
    fprintf(ff,'REMARK @     Ionic Strenght (M)                                         : %1.2f\n',IS);
    fprintf(ff,'REMARK @     Block Size                                                 : %d\n',bs);
    fprintf(ff,'REMARK %% Title "Details Regarding All Microstates Considered"\n');
    fprintf(ff,'REMARK %% Legend \n');
    fprintf(ff,'REMARK %%    column 1 "Number of Structured Blocks in the Microstate"  \n');
    fprintf(ff,'REMARK %%    column 2 "Statistical Weight of the Microstate"  \n');
    fprintf(ff,'REMARK %%    column 3 "Index of the First Block of the First Structured Region"  \n');
    fprintf(ff,'REMARK %%    column 4 "Length of the First Structured Region"  \n');
    fprintf(ff,'REMARK %%    column 5 "Index of the First Block of the Second Structured Region"  \n');
    fprintf(ff,'REMARK %%    column 6 "Length of the Second Structured Region"  \n');
    fprintf(ff,'REMARK %%    column 7 "Value Specifying the Approximation: 1 for SSA, 2 for DSA, 3 for DSAw/L"  \n');
    for i=1:size(pepval,1)
        fprintf(ff,'%3d %1.3e %3d %3d %3d %3d %1d\n',pepval(i,1),pepval(i,2),pepval(i,3),pepval(i,4),pepval(i,5),pepval(i,6),pepval(i,7));
    end
    fclose(ff);

    % Convert to pepval in desired format for further processing
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', [aa '_pepval.mat']);
    save(opPath, 'pepval');

    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'BlockApproximation.txt');
    ff = fopen(opPath,'wt');
    fprintf(ff,'REMARK @ Input Parameters\n');
    fprintf(ff,'REMARK @     van der Waals interaction energy (J/mol per native contact): %f\n',ene*1000);
    fprintf(ff,'REMARK @     Entropic Cost (J/mol.K per residue)                        : %3.2f\n',DS*1000);
    fprintf(ff,'REMARK @     Heat Capacity Change (J/mol.K per native contact)          : %3.2f\n',DCp*1000);
    fprintf(ff,'REMARK @     Temperature (K)                                            : %d\n',T);
    fprintf(ff,'REMARK @     Ionic Strenght (M)                                         : %1.2f\n',IS);
    fprintf(ff,'REMARK @     Block Size                                                 : %d\n',bs);
    fprintf(ff,'REMARK %% Title "Residue-to-Block Conversion"\n');
    fprintf(ff,'REMARK %% Legend \n');
    fprintf(ff,'REMARK %%    column 1 "Residue Index"  \n');
    fprintf(ff,'REMARK %%    column 2 "Block Index"  \n');
    for i=1:size(BlockDet,1)
        fprintf(ff,'%3d %3d\n',BlockDet(i,1),BlockDet(i,2));
    end
    fclose(ff);
    
    %% Generating images
    % 1D FE Profile
    fig1 = figure("Visible", false);
    set(fig1, 'Position', [1 1 1000 500]);
    plot(fes(:,1), fes(:,2), 'k', 'LineWidth', 2); axis([0 inf 0 inf]);
    if bs == 1
        xlabel('# of Structured Residues');
    else
        xlabel('# of Structured Blocks');
    end
    ylabel('Free Energy (kJ/mol)');
    xt = 0:10:size(fes,1); xt1 = {}; for i=1:length(xt); xt1{i}=num2str(xt(i)); end;
    set(gca, 'XTick', xt, 'XTickLabel', xt1);
    yt = 0:20:max(fes(:,2)); yt1 = {}; for i=1:length(yt); yt1{i}=num2str(yt(i)); end;
    set(gca, 'YTick', yt, 'YTickLabel', yt1);
    set(gca,'fontweight','bold','TickDir','out', 'LineWidth', 2, 'FontSize', 18)
    set(gca,'box','off','color','none')
    b = axes('Parent', fig1, 'Position',get(fig1,'Position'),'box','on','xtick',[],'ytick',[]);
    set(b,'LineWidth', 2)
    axes(fig1)
    linkprop([fig1 b],'position');
    xlim([fes(1,1), fes(end,1)])
    set(gcf,'PaperPositionMode','auto');
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'image1');
    print(char(opPath), '-dpng', '-r300')
    
    % 2D FE Profile
    fig2 = figure("Visible", false);
    set(gcf, 'Position', get(0, 'Screensize'));
    surf((0:nres-hv),(0:hv),fes2D); shading interp; colormap jet;
    axis tight
    xlabel('n_{C-term}');
    ylabel('n_{N-term}');
    zlabel('Free Energy (kJ/mol)');
    xt = 0:10:size(fes2D,2); xt1 = {}; for i=1:length(xt); xt1{i}=num2str(xt(i)); end;
    set(gca, 'XTick', xt, 'XTickLabel', xt1);
    yt = 0:10:size(fes2D,1); yt1 = {}; for i=1:length(yt); yt1{i}=num2str(yt(i)); end;
    set(gca, 'YTick', yt, 'YTickLabel', yt1);
    zt = 20:20:max(fes2D(~isinf(fes2D))+20); zt1 = {}; for i=1:length(zt); zt1{i}=num2str(zt(i)); end
    set(gca, 'ZTick', zt, 'ZTickLabel', zt1);
    pbaspect([1 1 0.5])
    grid off;
    set(gca,'fontweight','bold','TickDir','out', 'FontSize', 18)
    view([138 19])
    set(gcf,'PaperPositionMode','auto');
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'image2');
    print(char(opPath), '-dpng', '-r300')
    
    % Residue folding probability vs reaction coordinate
    fig3 = figure("Visible", false);
    set(gcf, 'Position', [1 1 1000 500]);
    pcolor(Fpath); colormap jet; shading interp;
    if bs == 1
        xlabel('# of Structured Residues');
    else
        xlabel('# of Structured Blocks');
    end
    ylabel('Residue Index');
    cb = colorbar; cb.Label.String = 'Probability'; %cb.Label.Rotation = -90; cb.Label.Position = cb.Label.Position+[2.25 0 0];
    xt = 0:10:size(Fpath,2); xt1 = {}; for i=1:length(xt); xt1{i}=num2str(xt(i)); end
    set(gca, 'XTick', xt, 'XTickLabel', xt1);
    if nres1 <= 200
        yspacing = 20;
    else
        yspacing = round(nres1/10,-1);
    end
    yt = 0:yspacing:size(Fpath,1); yt1 = {}; for i=1:length(yt); yt1{i}=num2str(yt(i)); end
    set(gca, 'YTick', yt, 'YTickLabel', yt1);
    set(gca,'fontweight','bold','TickDir','out', 'LineWidth', 2, 'FontSize', 18)
    %xl = get(gca,'XLabel'); xlp = get(xl, 'Position'); xlp(2) = xlp(2)+1; set(xl, 'Position', xlp);
    set(gca,'box','off','color','none')
    b = axes('Parent', fig3, 'Position',get(fig3,'Position'),'box','on','xtick',[],'ytick',[]);
    set(b,'LineWidth', 2)
    axes(fig3)
    linkprop([fig3 b],'position');
    set(gcf,'PaperPositionMode','auto');
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'image3');
    print(char(opPath), '-dpng', '-r300')
    
    % 2D FE Profile - Flat
    fig4 = figure("Visible", false);
    set(gcf, 'Position', [1 1 1000 1000]);
    pcolor((0:nres-hv),(0:hv),fes2D); shading interp; colormap jet;
    xlabel('n_{C-term}');
    ylabel('n_{N-term}');
    xt = 0:10:size(fes2D,2); xt1 = {}; for i=1:length(xt); xt1{i}=num2str(xt(i)); end
    set(gca, 'XTick', xt, 'XTickLabel', xt1);
    yt = 0:10:size(fes2D,1); yt1 = {}; for i=1:length(yt); yt1{i}=num2str(yt(i)); end
    set(gca, 'YTick', yt, 'YTickLabel', yt1);
    cb = colorbar; cb.Label.String = 'FE (kJ/mol)';
    set(gca,'fontweight','bold','TickDir','out', 'FontSize', 18)
    set(gcf,'PaperPositionMode','auto');
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'image4');
    print(char(opPath), '-dpng', '-r300')    

    % Number of states from each approximation
    fig5 = figure("Visible", false); hold on;
    set(gcf, 'Position', [1 1 600 500]);
    color_SSA = [141,160,203]/255;
    color_DSA = [252,141,98]/255;
    color_DSAwL = [102,194,165]/255;
    numstates_SSA = nnz(pepval(:,7)==1);
    numstates_DSA = nnz(pepval(:,7)==2);
    numstates_DSAwL = nnz(pepval(:,7)==3);
    bar(1, numstates_SSA, 'FaceColor', color_SSA);
    bar(2, numstates_DSA, 'FaceColor', color_DSA);
    bar(3, numstates_DSAwL, 'FaceColor', color_DSAwL);
    set(gca, 'YScale', 'log')
    ax = gca;
    ax.TickLength=[0.02 0.025];
    set(gca,'XTick',[1 2 3])
    set(gca,'xticklabel',{'SSA','DSA','DSAw/L'})
    ylabel('# of Microstates')
    set(gca,'fontweight','bold','FontSize', 18)
    set(gcf,'PaperPositionMode','auto');
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'image5');
    print(char(opPath), '-dpng', '-r300')
    
    % Partition function for each approximation
    fig6 = figure("Visible", false); hold on;
    set(gcf, 'Position', [1 1 600 500]);
    states_SSA = pepval(:,7)==1; pf_SSA = sum(pepval(states_SSA,2))/Zfin;
    states_DSA = pepval(:,7)==2; pf_DSA = sum(pepval(states_DSA,2))/Zfin;
    states_DSAwL = pepval(:,7)==3; pf_DSAwL = sum(pepval(states_DSAwL,2))/Zfin;
    approximations = categorical({'SSA','DSA','DSAw/L'});
    approximations = reordercats(approximations,{'SSA','DSA','DSAw/L'});
    bar(1, pf_SSA*100, 'FaceColor', color_SSA);
    bar(2, pf_DSA*100, 'FaceColor', color_DSA);
    bar(3, pf_DSAwL*100, 'FaceColor', color_DSAwL);
    ax = gca;
    ax.TickLength=[0.02 0.025];
    ylim([0 100])
    ylabel('% Contribution')
    set(gca,'XTick',[1 2 3])
    set(gca,'xticklabel',{'SSA','DSA','DSAw/L'})
    set(gca,'fontweight','bold','FontSize', 18)
    set(gcf,'PaperPositionMode','auto');
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'image6');
    print(char(opPath), '-dpng', '-r300')
    
    opPath = fullfile(pwd, '..', 'data', 'WSME_outputs', 'Statistics.txt');
    ff = fopen(opPath,'wt');
    fprintf(ff,'REMARK @ Input Parameters\n');
    fprintf(ff,'REMARK @     van der Waals interaction energy (J/mol per native contact): %f\n',ene*1000);
    fprintf(ff,'REMARK @     Entropic Cost (J/mol.K per residue)                        : %3.2f\n',DS*1000);
    fprintf(ff,'REMARK @     Heat Capacity Change (J/mol.K per native contact)          : %3.2f\n',DCp*1000);
    fprintf(ff,'REMARK @     Temperature (K)                                            : %d\n',T);
    fprintf(ff,'REMARK @     Ionic Strenght (M)                                         : %1.2f\n',IS);
    fprintf(ff,'REMARK @     Block Size                                                 : %d\n',bs);
    fprintf(ff,'REMARK %% Title "Number of Microstates from Each Approximation"\n');
    fprintf(ff,'SSA    : %d\n',numstates_SSA);
    fprintf(ff,'DSA    : %d\n',numstates_DSA);
    fprintf(ff,'DSAw/L : %d\n',numstates_DSAwL);
    fprintf(ff,'% Contribution from SSA to Total Partition Function  : %2.1f\n',pf_SSA*100);
    fprintf(ff,'% Contribution from DSA to Total Partition Function  : %2.1f\n',pf_DSA*100);
    fprintf(ff,'% Contribution from DSAw/L to Total Partition Function  : %2.1f\n',pf_DSAwL*100);
    fclose(ff);

    fclose('all');
end