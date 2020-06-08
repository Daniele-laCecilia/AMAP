function [AMAP, AMAE, AMAV, AMAskew, AMAkurt, Si, cond_mean, cond_var, cond_skew, cond_kurt ,cond_PThr, number, lower, upper] = gsa_amaP_DLC(sampling_points, Output_Mat, nclass,thr,name_var )
%

% this function computes AMA sensitivity indices and principal Sobol Indices (Si)
% starting from a set of Monte Carlo realizations. It also gives
% conditional moments with respect to one parameter at a time up to order 4
% inputs:
% sampling_points: MC realizations of the parameters [Nreal x Npar]
% Output_Mat: Full model evaluations [Nreal x NOutputs]
% nclass: nu,ber of classes used to compute conditional statistics
% thr: threshold value
% name_var: name of each output variable (text string)



ntime = size(Output_Mat,2);
Npar = size(sampling_points,2);
Nsim = size(sampling_points,1);

lower= zeros(Npar,nclass);
upper = zeros(Npar,nclass);
cond_mean = zeros(Npar,nclass);
cond_var = zeros(Npar,nclass);
cond_skew = zeros(Npar,nclass);
cond_kurt = zeros(Npar,nclass);

number = zeros(Npar,nclass);
cond_PThr = zeros(Npar,nclass);

AMAP = zeros(Npar,ntime);
AMAE = zeros(Npar,ntime);
AMAV = zeros(Npar,ntime);
Si = zeros(Npar,ntime);
AMAskew = zeros(Npar,ntime);
AMAkurt = zeros(Npar,ntime);

for itime = 1:ntime
    output = Output_Mat(:,itime);
    
    Pthr_unc = length(find(output>thr))/Nsim;
    
    for ipar = 1:Npar
        param_i = sampling_points(:,ipar);
        [sort_par,isort] = sort(param_i);
        for iclass = 1:nclass
            
            
            lower(ipar,iclass) = param_i(isort(round(Nsim*(iclass-1)/nclass)+1));
            upper(ipar,iclass) = param_i(isort(round(Nsim*(iclass)/nclass)));
            
            j_sampling_points_class = find(param_i>lower(ipar,iclass) & param_i<upper(ipar,iclass));
            output_class = output(j_sampling_points_class);
            cond_mean(ipar,iclass) = mean(output_class);
            cond_var(ipar,iclass) = var(output_class);
            cond_skew(ipar,iclass) = skewness(output_class);
            cond_kurt(ipar,iclass) = kurtosis(output_class);
            
            number(ipar,iclass) = length(j_sampling_points_class);
            cond_PThr(ipar,iclass) = length(find(output_class>thr))/number(ipar,iclass) ;
        end
        
        AMAP(ipar,itime) = mean(abs(cond_PThr(ipar,:)-Pthr_unc));
        AMAE(ipar,itime) = mean(abs(cond_mean(ipar,:)-mean(output)))/abs(mean(output));
        AMAV(ipar,itime) = mean(abs(cond_var(ipar,:)-var(output)))/var(output);
        Si(ipar,itime) =  var(cond_mean(ipar,:))/var(output);
        AMAskew(ipar,itime) = mean(abs(cond_skew(ipar,:)-skewness(output)))/abs(skewness(output));
        AMAkurt(ipar,itime) = mean(abs(cond_kurt(ipar,:)-kurtosis(output)))/kurtosis(output);
    end
end

end
