function diffuse_fraction_test(MD,SP,dt)
% DIFFUSE_FRACTION_TEST(MD,SUNEL,ITS) - Test all diffuse-fraction models in DIFFUSE_FRACTION on
%   meteo-data MD, and generate a box-plot of model errors.
%
% INPUT:
%   MD structure with equal-sized vector fields {GHI [W/m²], DHI [W/m²], kt, Patm [Pa], and
%       tpw [kg/m²]}, see COMPLETEMETEODATA for details.
%   SUNEL vector of solar elevation angles (degrees) for each data point.
%   ITS - boolean scalar, passed directly to DIFFUSE_FRACTION.
%
% See also: DIFFUSE_FRACTION, COMPLETEMETEODATA

    MDL = {'Orgill','Erbs','Reindl1','Reindl2','DISC','SOT2','SOT3','SOT4','DIRINT','Engerer2'};
    n = numel(MDL);

    ITS = nargin > 3 && ~isempty(dt);
    parsestruct(MD,{'GHI','DHI','Patm','tpw','CSGHI','ENI'},'-n','-v','-r','-e');
    parsestruct(SP,{'El'},'opt',{'w'},'-n','-v','-r','size',size(MD.GHI));

    notdark = SP.El > 0 & MD.GHI > 0;
    EHI = max(0,MD.ENI.*sind(SP.El));
    MD.kt = MD.GHI./EHI;
    MD.ktc = MD.CSGHI./EHI;
    
    warning_disabler = naptime('diffuse_fraction:ignored'); %#ok<NASGU>

    [Kd,Kn] = cellfun(@(mdl) diffuse_fraction(MD.kt(notdark),90-SP.El(notdark),MD.Patm(notdark),...
            MD.tpw(notdark),ITS,MD.ktc(notdark),SP.w(notdark),[],dt,mdl),MDL,'unif',0);
        
    MDL{end+1} = 'Ensemble'; n = n+1;
    w = 1./mean((MD.ENI(notdark).*[Kn{:}] - MD.BNI(notdark)).^2,1,'omitnan');
    w = w'./sum(w);
    Kd{end+1} = [Kd{:}]*w;
    Kn{end+1} = [Kn{:}]*w;
        
    err = cat(3, MD.GHI(notdark).*[Kd{:}] - MD.DHI(notdark),...
                 MD.ENI(notdark).*[Kn{:}] - MD.BNI(notdark));

    MBE = shiftdim(mean(err,1,'omitnan'),1);
    RMS = hypot(MBE,shiftdim(std(err,1,'omitnan'),1));
    MAD = shiftdim(mad(err,1,1),1);

    tags = arrayfun(@(a,b,c) sprintf('\nMBE:%0.1f\nRMS:%0.1f\nMAD:%0.1f',a,b,c),...
        MBE,RMS,MAD,'unif',0);

    GUIfigure('diffuse_fraction_test'); clf();
    for j = 1:2
        subplot(2,1,j);
        boxplot(err(:,:,j),'Labels',MDL);
        switch j
            case 1, ylabel('DHI_{model} - DHI_{meas} [W/m^2]');
            case 2, ylabel('DNI_{model} - DNI_{meas} [W/m^2]');
        end
        grid on;

        ymax = max(abs(ylim()));
        ylim([-1,1]*ymax);

        text(1:n,repmat(ymax,1,n),tags(:,j),...
                'fontsize',10,'verticalalignment','top','horizontalalignment','center');
    end
end
