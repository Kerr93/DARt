library(EnvStats)
#library(lubridate)
#library('spatialfil')
library(ggplot2)

# main function
est_cutoff <- function(cutoff_range,sigma_range,theta_prior0,N2,pt_range,dt_range,beta0,It,R0,SI){
    res_all = vector(mode = "list", length = length(cutoff_range))
    t_start1 = 2
    t_end2 = N2
    for (n_c in 1:length(cutoff_range)){
        t_end1 = cutoff_range[n_c]-1
        t_start2 = cutoff_range[n_c]
        res_piece1 = piece_estimate(t_start1,t_end1,sigma_range,theta_prior0,N2,pt_range,dt_range,beta0,It)
        post_A_matrix1 = res_piece1$post_A_matrix
        prior_A_matrix1 = res_piece1$prior_A_matrix
        res_piece2 = piece_estimate(t_start2,t_end2,sigma_range,theta_prior0,N2,pt_range,dt_range,beta0,It)
        post_A_matrix2 = res_piece2$post_A_matrix
        prior_A_matrix2 = res_piece2$prior_A_matrix

        res_all[[n_c]] = list("post_A_matrix1"=  post_A_matrix1,
                              "post_A_matrix2" = post_A_matrix2,
                              "prior_A_matrix1" = prior_A_matrix1,
                              "prior_A_matrix2" = prior_A_matrix2)

    }

    PROB_H <- cal_high_level_evidence(cutoff_range,sigma_range,res_all)
    PROB_H_norm = PROB_H/sum(PROB_H)
  
    # for low level distributions
    res_low <- cal_low_level_evidence(PROB_H_norm,cutoff_range,sigma_range,res_all,N2,pt_range,dt_range)
    p_d_prior = res_low$p_d_prior
    p_d_prob = res_low$p_d_prob
    
    Rt_samples <- Rt_samples_gen(pt_range,dt_range,SI,N2,R0)
 
    res <- summarize_res(p_d_prob,p_d_prior,pt_range,dt_range,N2,Rt_samples,It,PROB_H)
    return(res)
}




summarize_res <- function(p_d_prob,p_d_prior,pt_range,dt_range,N2,Rt_samples, It,PROB_H){
    date = 1:N2
    p_t_mean = rep(0,N2)
    d_t_mean = rep(0,N2)
    I_t_mean = rep(0,N2)
    R_t_mean = rep(0,N2)
    I_pred_t_mean = rep(0,N2)
     # It_samples will be used for calculating I
    It_samples = array(0,dim=c(N2,length(pt_range),length(dt_range)))
    for (i in 2:N2){
        I_i = cal_I_sample(pt_range,dt_range,beta0,It,i,N2)
        It_samples[i,,] = I_i$I_pred_samples
    }
    
    for (i in 2:N2){
        pd_i = p_d_prob[i,,]/sum(p_d_prob[i,,])
        pd_prior_i = p_d_prior[i,,]/sum(p_d_prior[i,,])
        for(n_p in 1:length(pt_range)){
            p  = pt_range[n_p]
            for(n_d in 1:length(dt_range)){
                d = dt_range[n_d]
                p_t_mean[i] = p_t_mean[i]+p*pd_i[n_p,n_d]
                d_t_mean[i] = d_t_mean[i]+d*pd_i[n_p,n_d]
                I_t_mean[i] = I_t_mean[i]+It_samples[i,n_p,n_d]*pd_i[n_p,n_d]
                R_t_mean[i] = R_t_mean[i]+Rt_samples[n_p,n_d]*pd_i[n_p,n_d]
                I_pred_t_mean[i] = I_pred_t_mean[i] + It_samples[i,n_p,n_d]*pd_prior_i[n_p,n_d]

            }
        }
    } 
    pt25 = rep(0,N2)
    pt75 = rep(0,N2)
    dt25 = rep(0,N2)
    dt75 = rep(0,N2)

    Rt25 = rep(0,N2)
    Rt75 = rep(0,N2)
    for (i in 2:N2){
        p_i = apply(p_d_prob[i,,], 1, sum)
        pt25[i] = pt_range[which(cumsum(p_i)>0.025)[1]]
        pt75[i] = pt_range[which(cumsum(p_i)>=0.975)[1]]

        d_i = apply(p_d_prob[i,,], 2, sum)
        dt25[i] = dt_range[which(cumsum(d_i)>0.025)[1]]
        dt75[i] = dt_range[which(cumsum(d_i)>=0.975)[1]]

        xi = sort(as.vector(Rt_samples), index.return=TRUE)
        R_val = xi$x
        R_oder = xi$ix
        R_prob = cumsum(as.vector(p_d_prob[i,,])[R_oder])
        Rt25[i] = R_val[which(R_prob>0.025)[1]]
        Rt75[i] = R_val[which(R_prob>=0.975)[1]]
    }

    result_epidemic = data.frame('date' = date, 'It'= It, 'I_correct_all' = I_t_mean, 'I_pred_all' = I_pred_t_mean, 
                                 'R_t_mean' = R_t_mean,'Rt025' = Rt25, 'Rt975'= Rt75,
                                 'p_t_mean' = p_t_mean, 'pt025'= pt25, 'pt975'= pt75,
                                 'd_t_mean' = d_t_mean, 'dt025'= dt25, 'dt975'= dt75)

    # distribution of cutoff date
    cutoff_vec = apply(PROB_H, 1, sum)
    cutoff_vec_prob = cutoff_vec/norm(cutoff_vec)
#     which.max(cutoff_vec_prob)
#     date[cutoff_range[which.max(cutoff_vec_prob)]]
    cutoff_dist = data.frame('date'=date[cutoff_range],'prob' = cutoff_vec_prob)


    # distribution of two sigma values

    sigma1 = apply(PROB_H, 2, sum)
    sigma1_norm = sigma1/sum(sigma1)

    sigma2 = apply(PROB_H, 3, sum)
    sigma2_norm = sigma2/sum(sigma2)

    sigma12_marginal_dist = data.frame('Sigma'=sigma_range,'prob_sigma1' = sigma1_norm,'prob_sigma2' = sigma2_norm)


    sigma12 = apply(PROB_H, c(2,3), sum)
    sigma12_norm = sigma12/sum(sigma12)
    colnames(sigma12_norm) = sigma_range
    rownames(sigma12_norm) = sigma_range
    
    return(list('result_epidemic'=result_epidemic, 'cutoff_dist'=cutoff_dist,
                'sigma12_marginal_dist'= sigma12_marginal_dist, 'sigma12_norm' = sigma12_norm))
    
}

serial_interval <-  function(mean1,cv1,N2){
    # discritized serial interval
    x1 = rgammaAlt(5e6,mean1,cv1)
    f = ecdf(x1)
    h = rep(0,N2)
    # h = rep(0,forecast+N)
    h[1] = (f(1.5) - f(0)) 
    for(i in 2:length(h)) {
      h[i] = (f(i+.5) - f(i-.5)) / (1-f(i-.5))
    }
    s = rep(0,N2)
    s[1] = 1 
    for(i in 2:N2) {
        s[i] = s[i-1]*(1-h[i-1])
    }
    SI0 = s * h  
    SI1 = SI0[!is.nan(SI0)]
    SI = rep(0,N2)
    SI[1:length(SI1)] = SI1
    
    return(SI)
}

discrete_weibull <-  function(shape,scale,N2){
# https://rdrr.io/r/stats/Weibull.html
    x1 = rweibull(5e6, shape, scale)
    f = ecdf(x1)
    h = rep(0,N2)
    # h = rep(0,forecast+N)
    h[1] = (f(1.5) - f(0)) 
    for(i in 2:length(h)) {
      h[i] = (f(i+.5) - f(i-.5)) / (1-f(i-.5))
    }
    s = rep(0,N2)
    s[1] = 1 
    for(i in 2:N2) {
        s[i] = s[i-1]*(1-h[i-1])
    }
    SI0 = s * h  
    SI1 = SI0[!is.nan(SI0)]
    SI = rep(0,N2)
    SI[1:length(SI1)] = SI1
    
    return(SI)
}

discrete_LogNormal <-  function(logmean,logsd,N2){
    # https://rdrr.io/cran/EnvStats/man/elnorm.html
    mean = exp(logmean+logsd^2/2)
    cv = sqrt(exp(logsd^2)-1)
    x1 = rlnormAlt(5e6, mean, cv )
    f = ecdf(x1)
    h = rep(0,N2)
    # h = rep(0,forecast+N)
    h[1] = (f(1.5) - f(0)) 
    for(i in 2:length(h)) {
      h[i] = (f(i+.5) - f(i-.5)) / (1-f(i-.5))
    }
    s = rep(0,N2)
    s[1] = 1 
    for(i in 2:N2) {
        s[i] = s[i-1]*(1-h[i-1])
    }
    SI0 = s * h  
    SI1 = SI0[!is.nan(SI0)]
    SI = rep(0,N2)
    SI[1:length(SI1)] = SI1
    
    return(SI)
}
discrete_Gamma <-  function(mean1,cv1,N2){

#     https://www.rdocumentation.org/packages/EnvStats/versions/2.3.1/topics/GammaAlt
#     https://rdrr.io/cran/extraDistr/man/DiscreteGamma.html
    x1 = rgammaAlt(5e6,mean1,cv1)
    f = ecdf(x1)
    h = rep(0,N2)
    # h = rep(0,forecast+N)
    h[1] = (f(1.5) - f(0)) 
    for(i in 2:length(h)) {
      h[i] = (f(i+.5) - f(i-.5)) / (1-f(i-.5))
    }
    s = rep(0,N2)
    s[1] = 1 
    for(i in 2:N2) {
        s[i] = s[i-1]*(1-h[i-1])
    }
    SI0 = s * h  
    SI1 = SI0[!is.nan(SI0)]
    SI = rep(0,N2)
    SI[1:length(SI1)] = SI1
    
    return(SI)
}


T_A = function(u,sigma) (applyFilter(x = u, kernel = convKernel(sigma = sigma, k = 'gaussian'))) # identical to prior (need to 
    
# cal_low_level_evidence <- function(PROB_H_norm,cutoff_range,sigma_range,res_all,N2,pt_range,dt_range){

#     p_d_prob = array(0,dim=c(N2,length(pt_range),length(dt_range)))
#     p_d_prior = array(0,dim=c(N2,length(pt_range),length(dt_range)))
#     for (n_c in 1:length(cutoff_range)){
#         t_end1 = cutoff_range[n_c]-1
#         post_A_matrix1 = res_all[[n_c]]$post_A_matrix1
#         post_A_matrix2 = res_all[[n_c]]$post_A_matrix2

#         prior_A_matrix1 = res_all[[n_c]]$prior_A_matrix1
#         prior_A_matrix2 = res_all[[n_c]]$prior_A_matrix2

#         for(n_s1 in 1:length(sigma_range)){
#                 for(n_s2 in 1:length(sigma_range)){
#                     probH = PROB_H_norm[n_c,n_s1,n_s2]
#                     cutoff_prior = post_A_matrix1[n_s1,,,t_end1]/sum(post_A_matrix1[n_s1,,,t_end1])
#                     for (i in 2:N2){                   
#                         if (i <= t_end1){
#                             post1_norm_i = post_A_matrix1[n_s1,,,i]/sum(post_A_matrix1[n_s1,,,i])
#                             prior1_norm_i = prior_A_matrix1[n_s1,,,i]/sum(prior_A_matrix1[n_s1,,,i])
#                             for(n_p in 1:length(pt_range)){
#                                 for(n_d in 1:length(dt_range)){
#                                      p_d_prob[i,n_p,n_d] = p_d_prob[i,n_p,n_d] + post1_norm_i[n_p,n_d]*probH

#                                      p_d_prior[i,n_p,n_d] = p_d_prior[i,n_p,n_d] + prior1_norm_i[n_p,n_d]*probH
#                                 }
#                             }
#                         }else{
#                             post2_norm_i = post_A_matrix2[n_s2,,,i]/sum(post_A_matrix2[n_s2,,,i])
#                             prior2_norm_i = prior_A_matrix2[n_s2,,,i]/sum(prior_A_matrix2[n_s2,,,i])

#                             for(n_p in 1:length(pt_range)){
#                                 for(n_d in 1:length(dt_range)){
#                                     p_d_prob[i,n_p,n_d] = p_d_prob[i,n_p,n_d] + post2_norm_i[n_p,n_d]*probH
#                                     if (i == t_end1+1){
#                                         p_d_prior[i,n_p,n_d] = p_d_prior[i,n_p,n_d] + cutoff_prior[n_p,n_d]*probH
#                                     }else{
#                                         p_d_prior[i,n_p,n_d] = p_d_prior[i,n_p,n_d] + prior2_norm_i[n_p,n_d]*probH
#                                     }

#                                 }
#                             }                
#                         }


#                     }

#                 }
#             }    


#       }
#     return(list('p_d_prior' = p_d_prior, 'p_d_prob' = p_d_prob))
# }
cal_low_level_evidence <- function(PROB_H_norm,cutoff_range,sigma_range,res_all,N2,pt_range,dt_range){

    p_d_prob = array(0,dim=c(N2,length(pt_range),length(dt_range)))
    p_d_prior = array(0,dim=c(N2,length(pt_range),length(dt_range)))
    for (n_c in 1:length(cutoff_range)){
        t_end1 = cutoff_range[n_c]-1
        post_A_matrix1 = res_all[[n_c]]$post_A_matrix1
        post_A_matrix2 = res_all[[n_c]]$post_A_matrix2

        prior_A_matrix1 = res_all[[n_c]]$prior_A_matrix1
        prior_A_matrix2 = res_all[[n_c]]$prior_A_matrix2

        for(n_s1 in 1:length(sigma_range)){
                for(n_s2 in 1:length(sigma_range)){
                    probH = PROB_H_norm[n_c,n_s1,n_s2]
                    if (sum(post_A_matrix1[n_s1,,,t_end1])==0){
                        cutoff_prior = post_A_matrix1[n_s1,,,t_end1]
                    }else{
                        cutoff_prior = post_A_matrix1[n_s1,,,t_end1]/sum(post_A_matrix1[n_s1,,,t_end1])
                    }
                    
                    for (i in 2:N2){                   
                        if (i <= t_end1){
                            if (sum(post_A_matrix1[n_s1,,,i])==0){
                                post1_norm_i = post_A_matrix1[n_s1,,,i]
                            }else{
                                post1_norm_i = post_A_matrix1[n_s1,,,i]/sum(post_A_matrix1[n_s1,,,i])
                            }
                            if (sum(prior_A_matrix1[n_s1,,,i])==0){
                                prior1_norm_i = prior_A_matrix1[n_s1,,,i]
                            }else{
                                prior1_norm_i = prior_A_matrix1[n_s1,,,i]/sum(prior_A_matrix1[n_s1,,,i])
                            }                            
                            
                            for(n_p in 1:length(pt_range)){
                                for(n_d in 1:length(dt_range)){
                                     p_d_prob[i,n_p,n_d] = p_d_prob[i,n_p,n_d] + post1_norm_i[n_p,n_d]*probH

                                     p_d_prior[i,n_p,n_d] = p_d_prior[i,n_p,n_d] + prior1_norm_i[n_p,n_d]*probH
                                }
                            }
                        }else{
                            if (sum(post_A_matrix2[n_s2,,,i])==0){
                                post2_norm_i = post_A_matrix2[n_s2,,,i]
                            }else{
                                post2_norm_i = post_A_matrix2[n_s2,,,i]/sum(post_A_matrix2[n_s2,,,i])
                            }
                            if (sum(prior_A_matrix2[n_s2,,,i])==0){
                                prior2_norm_i = prior_A_matrix2[n_s2,,,i]
                            }else{
                                prior2_norm_i = prior_A_matrix2[n_s2,,,i]/sum(prior_A_matrix2[n_s2,,,i])
                            } 
                          
                            for(n_p in 1:length(pt_range)){
                                for(n_d in 1:length(dt_range)){
                                    p_d_prob[i,n_p,n_d] = p_d_prob[i,n_p,n_d] + post2_norm_i[n_p,n_d]*probH
                                    if (i == t_end1+1){
                                        p_d_prior[i,n_p,n_d] = p_d_prior[i,n_p,n_d] + cutoff_prior[n_p,n_d]*probH
                                    }else{
                                        p_d_prior[i,n_p,n_d] = p_d_prior[i,n_p,n_d] + prior2_norm_i[n_p,n_d]*probH
                                    }

                                }
                            }                
                        }


                    }

                }
            }    


      }
    return(list('p_d_prior' = p_d_prior, 'p_d_prob' = p_d_prob))
}



cal_high_level_evidence <- function(cutoff_range,sigma_range,res_all){
    # get high level distribution
    PROB_H = array(0,dim=c(length(cutoff_range),length(sigma_range),length(sigma_range))) 
    for (n_c in 1:length(cutoff_range)){
        t_end1 = cutoff_range[n_c]-1
        post_A_matrix1 = res_all[[n_c]]$post_A_matrix1
        post_A_matrix2 = res_all[[n_c]]$post_A_matrix2

        for(n_s1 in 1:length(sigma_range)){
                for(n_s2 in 1:length(sigma_range)){
                    PROB_H[n_c,n_s1,n_s2] = mean(post_A_matrix1[n_s1,,,t_end1])*mean(post_A_matrix2[n_s2,,,N2])
                    #####################
#                     PROB_H[n_c,n_s1,n_s2] = 10^100*mean(post_A_matrix1[n_s1,,,t_end1])*mean(post_A_matrix2[n_s2,,,N2])
                    #####################
                }
            }    


        }
    
    return(PROB_H)
}



piece_estimate <-  function(t_start, t_end, sigma_range,theta_prior0,N2,pt_range,dt_range,beta0,It){
 
    post_A_matrix = array(0,dim=c(length(sigma_range),length(pt_range),length(dt_range),N2))  
    prior_A_matrix = array(0,dim=c(length(sigma_range),length(pt_range),length(dt_range),N2)) 
    
    for(n_s in 1:length(sigma_range)){
        post_theta = theta_prior0
        sigma = sigma_range[n_s]
    
        for(i in t_start:t_end){         
            # generate I samples; to be used for likelihood calculation
            res = cal_I_sample(pt_range,dt_range,beta0,It,i,N2)
            prob = res$prob
     
            # prior
            pior_theta_A = T_A(post_theta,sigma)# smooth change

            # unnormalized posterior for model evidence calculation
            ProbA = 10^2*prob*pior_theta_A##########5
#             # normalized posterior for next step prior definition
#             post_theta = ProbA/sum(ProbA) 
            post_theta = ProbA
            # save evidence
            post_A_matrix[n_s,,,i] = ProbA
            
            prior_A_matrix[n_s,,,i] = pior_theta_A##########
        }

    }
    return(list("post_A_matrix" = post_A_matrix, "prior_A_matrix" = prior_A_matrix))
}



Rt_samples_gen <-  function(pt_range,dt_range,SI,N2,R0){
    # Rt_samples will be used for calculating Rt
    Rt_samples = matrix(0,length(pt_range), length(dt_range))
    for(n_p in 1:length(pt_range)){
        for(n_d in 1:length(dt_range)){
            p = pt_range[n_p]
            d = dt_range[n_d]
            SIt = SI
            SIt[seq(1,N2,1)> d] = 0 
            Rt_samples[n_p,n_d] = R0*p*sum(SIt)
        }
    }
    return(Rt_samples)
}
cal_I_sample <-  function(pt_range,dt_range,beta0,It,i,N2){
    #likelihood of I(t): I_pred_samples & prob
    # beta_t distribution
    beta_t_samples = matrix(0, length(pt_range)*length(dt_range), length(beta0))
    n_sample = 0
    I_pred_samples = matrix(0, length(pt_range),length(dt_range))
    prob = matrix(0, length(pt_range),length(dt_range))
    for(n_p in 1:length(pt_range)){
        for(n_d in 1:length(dt_range)){
            n_sample = n_sample+1
            p = pt_range[n_p]
            d = dt_range[n_d]
            beta_t = beta0*p
            beta_t[seq(1,N2,1)> d] = 0 ############need to check >= or >
            beta_t_samples[n_sample,] = beta_t
            I_pred = 0
            for(tao in 1:(i-1)){
                I_pred = I_pred + beta_t[i-tao]*It[tao]#beta_t???
            }
            
            I_pred_samples[n_p,n_d] = I_pred
            prob[n_p,n_d] = dpois(round(It[i]),I_pred)# record 'prob'
        }
    }
    res <- list("I_pred_samples" = I_pred_samples, "prob" = prob)
    return(res)
}

save_results <- function(city_name, res4){
    result_epidemic = res4$result_epidemic
    cutoff_dist = res4$cutoff_dist
    sigma12_marginal_dist = res4$sigma12_marginal_dist
    sigma12_norm = res4$sigma12_norm
    
    dir.create(paste('result_us/',city_name,sep=""))
    result_epidemic$date = date
    cutoff_dist$date = date[cutoff_range]
    write.csv(result_epidemic,paste('result_us/',city_name,'/',city_name,'_result_epidemic.csv',sep=""))
    write.csv(cutoff_dist,paste('result_us/',city_name,'/',city_name,'_cutoff_dist.csv',sep=""))
    write.csv(sigma12_marginal_dist,paste('result_us/',city_name,'/',city_name,'_sigma12_marginal_dist.csv',sep=""))
    

    pdf(paste('result_us/',city_name,'/',city_name,'_plot.pdf',sep=""))
    plot1 <- ggplot(result_epidemic, aes(x=date, y=p_t_mean)) +
       geom_line(size=1, alpha=0.8) +
       geom_ribbon(aes(ymin=pt025, ymax=pt975), fill="blue", alpha=0.2) + labs(y = "p(t)")+
       ggtitle(city_name)

    plot2 <- ggplot(result_epidemic, aes(x=date, y=d_t_mean)) +
       geom_line(size=1, alpha=0.8) +
       geom_ribbon(aes(ymin=dt025, ymax=dt975), fill="blue", alpha=0.2) + labs(y = "d(t)")

    plot3 <- ggplot(result_epidemic, aes(x=date, y=R_t_mean)) +
       geom_line(size=1, alpha=0.8) +
       geom_ribbon(aes(ymin=Rt025, ymax=Rt975), fill="blue", alpha=0.2) + labs(y = "R(t)")

    I_all = data.frame(date,result_epidemic$I_correct_al,result_epidemic$I_pred_all,result_epidemic$It)
    Molten <- melt(I_all, id.vars = "date")
    plot4 <- ggplot(Molten, aes(x = date, y = value, colour = variable)) + geom_line()+scale_color_discrete(name = "It", labels = c("It estimated", "It predicted", "It"))

    plot5 <- ggplot(cutoff_dist, aes(x = date, y = prob)) +geom_bar(stat="identity", fill="steelblue")

    Molten <- melt(sigma12_marginal_dist, id.vars = "Sigma")
    plot6 <- ggplot(Molten, aes(x = Sigma, y = value, colour = variable)) + geom_line()+scale_color_discrete(name = "Sigma", labels = c("Sigma1", "Sigma2"))

    multiplot(plot1,plot2,plot3,plot4,plot5,plot6, cols=2)
    dev.off() 
}
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

calgt <- function(shape, scale){
    discrete_weibull(shape, scale, 20)
    }

cut <- function(threshold,dis){
    index = which(dis > threshold)
    dis0 = numeric(index[1] - 1)
    c(dis0, dis[index] / sum(dis[index]))
    }

calinc <- function(shape, scale){
    disinc = discrete_LogNormal(shape, scale, 20)
    #indexinc < - which(disinc > 0.1)
    #inc0 < - numeric(indexinc[1] - 1)
    #c(inc0, disinc[indexinc] / sum(disinc[indexinc]))
    }


calrep <- function(repmean, repsd, incmean, incsd){
    N2 = 30
    x1 = rgammaAlt(5e6, repmean, repsd)

    mean = exp(incmean + incsd ^ 2 / 2)
    cv = sqrt(exp(incsd ^ 2) - 1)

    x2 = rlnormAlt(5e6, mean, cv)

    f = ecdf(x1 + x2)
    h = rep(0, N2)
    h[1] = (f(1.5) - f(0))
    for (i in 2:length(h)) {
        h[i] = (f(i+.5) - f(i-.5)) / (1-f(i-.5))
    }
    s = rep(0, N2)
    s[1] = 1
    for (i in 2:N2) {
        s[i] = s[i - 1] * (1 - h[i - 1])
    }
    SI0 = s * h
    SI1 = SI0[! is.nan(SI0)]
    SI = rep(0, N2)
    SI[1:length(SI1)] = SI1

    indexSI <- which(SI>0.09)
    SI0 <- numeric(indexSI[1]-1)
    c(SI0,SI[indexSI])

    }
