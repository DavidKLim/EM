List score_info_g(int N, double ph, arma::mat mu, arma::vec y, arma::mat wts){
    double lambda = 1e-25;
    double score1 = 0, info1 = 0;
    double inv_ph = 1/ph + lambda;  // try adding lambda to stabilize (when phi very large)?
    ph=ph/(1+ph*1e-25);
    double muic, yi, wtsic, scoreic, infoic, inv_phMuic;
    
    int n=y.size();
    int k = wts.ncol();
    
    for(int i=0; i<n; i++){
        for(int c=0; c<k; c++){
            yi = y(i);
            muic = mu(i,c);
            wtsic = wts(i,c);
        
            inv_phMuic = inv_ph + muic;
        
            scoreic = wtsic * (R::digamma(inv_ph+yi) - R::digamma(inv_ph) + log(inv_ph) + 1 - log(inv_phMuic) - (yi+inv_ph)/(inv_phMuic));
            infoic = wtsic * (R::trigamma(inv_ph) + 2/(inv_phMuic) - R::trigamma(inv_ph+yi) - ph - (yi+inv_ph)/pow(inv_phMuic,2));
        
            score1 += scoreic;
            info1 += infoic;
        }
    }
    
    double score = score1 * (-inv_ph*inv_ph) + 2*lambda*ph;
    double info = info1 * pow(inv_ph,4) + 2*lambda;
    
    return List::create(score,info);
}


double phi_ml_g(arma::vec y, arma::mat mu, arma::mat wts, int limit, int trace){
    double eps = 0.0001220703; /* from R */
    double p0 = 0;
    double N = accu(wts);
    int n = y.size();
    int k = wts.ncol();
    arma::mat wtd_y(n,k);
    
    int equal=0;
    if(n==1){
        return(p0=0);
    } /*else {
        for(int i=0; i<n; i++){
            for(int c=0; c<k; c++){
                wtd_y(i,c)=y(i)*wts(i,c);
            }
        }
    }
    
    double sum_wtd_y = accu(wtd_y);
    for(int i=0; i<n; i++){
        for(int c=0; c<k, c++){
            if(wtd_y(i,c)==sum_wtd_y){
                return(p0=0);
            }
            if(i<n-1){
                if(wtd_y(i,c) == wtd_y(i+1,c)){
                    equal++;
                }
            }
        }
    }
    
    if(equal==(n-1)){
        return(p0=0);
    }  */ /* I DONT THINK THIS IS NECESSARY?? */
    if(accu(mu)<1e-13){
        return(p0=0);
    }
    
    for(int i=0; i<n; i++){
        for(int c=0; c<k; c++){
            double wtsic=wts(i,c), yi=y(i), muic=mu(i,c);
            p0 += wtsic * pow(yi/muic-1,2);
        }
    }
    p0 = p0/N;
    
    int it=0;
    double del=1;
    
    if(trace==1){
        Rprintf("phi_ml: iter %d 'phi = %f' \n",it,p0);
    }
    while(it < limit && fabs(del) > eps){
        it += 1;
        p0 = fabs(p0);
        List scoreinfo = score_info(N,p0,mu,y,wts);
        double score=scoreinfo[0], info=scoreinfo[1];
        del = score/info;
        p0 += del;
        if(trace==1){
            Rprintf("score: %f\n",score);
            Rprintf("info: %f\n",info);
            Rprintf("phi_ml: iter %d 'phi = %f' \n",it,p0);
        }
    }
    
    if(p0 < 0){
        p0 = 0;
        if(trace==1){
            Rprintf("estimate truncated at zero \n");
        }
    }
    if(it == limit && trace==1){
        Rprintf("iteration limit reached \n");
    }
    
    return(p0);
}
