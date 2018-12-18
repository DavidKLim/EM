//
//  M.cpp
//  
//
//  Created by KYUNG T LIM on 1/19/18.
//
//  * NOTE ALL INDICES START FROM 0 IN C++, NOT 1 *

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

List score_info_g(int N, double ph, arma::mat mu, arma::vec y, arma::mat wts){
    double lambda = 1e-25;
    double score1 = 0, info1 = 0;
    double inv_ph = 1/ph + lambda;  // try adding lambda to stabilize (when phi very large)?
    ph=ph/(1+ph*1e-25);
    double muic, yi, wtsci, scoreic, infoic, inv_phMuic;
    
    int n=y.size();
    int k = mu.n_cols;
    
    for(int i=0; i<n; i++){
        for(int c=0; c<k; c++){
            yi = y(i);
            muic = mu(i,c);
            wtsci = wts(c,i);
            
            inv_phMuic = inv_ph + muic;
            
            scoreic = wtsci * (R::digamma(inv_ph+yi) - R::digamma(inv_ph) + log(inv_ph) + 1 - log(inv_phMuic) - (yi+inv_ph)/(inv_phMuic));
            infoic = wtsci * (R::trigamma(inv_ph) + 2/(inv_phMuic) - R::trigamma(inv_ph+yi) - ph - (yi+inv_ph)/pow(inv_phMuic,2));
            
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
    int k = mu.n_cols;
    arma::mat wtd_y(n,k);
    
    /*int equal=0;*/
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
            double wtsci=wts(c,i), yi=y(i), muic=mu(i,c);
            p0 += wtsci * pow(yi/muic-1,2);
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
        List scoreinfo = score_info_g(N,p0,mu,y,wts);
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


List score_info(int N, double ph, arma::vec mu, arma::vec y, arma::rowvec wts){
    double lambda = 1e-25;
    double score1 = 0, info1 = 0;
    double inv_ph = 1/ph + lambda;  // try adding lambda to stabilize (when phi very large)?
    ph=ph/(1+ph*1e-25);
    double mui, yi, wtsi, scorei, infoi, inv_phMui;
    
    int n=y.size();
    
    for(int i=0; i<n; i++){
        yi = y(i);
        mui = mu(i);
        wtsi = wts(i);
        
        inv_phMui = inv_ph + mui;
        
        scorei = wtsi * (R::digamma(inv_ph+yi) - R::digamma(inv_ph) + log(inv_ph) + 1 - log(inv_phMui) - (yi+inv_ph)/(inv_phMui));
        infoi = wtsi * (R::trigamma(inv_ph) + 2/(inv_phMui) - R::trigamma(inv_ph+yi) - ph - (yi+inv_ph)/pow(inv_phMui,2));
        
        score1 += scorei;
        info1 += infoi;
        //Rprintf("Score samp %d: %f \n",i+1,scorei);
        //Rprintf("Info samp %d: %f \n",i+1,infoi);
        
        /*Rprintf("Samp %d trigamma(inv_ph): %f\n",i+1,R::trigamma(inv_ph));
        Rprintf("Samp %d 2/(inv_phMui): %f\n",i+1,2/(inv_phMui));
        Rprintf("Samp %d trigamma(inv_ph+yi): %f\n",i+1,R::trigamma(inv_ph+yi));
        Rprintf("Samp %d ph: %f\n",i+1,ph);
        Rprintf("Samp %d (yi+inv_ph)/inv_phMui^2: %f\n",i+1,(yi+inv_ph)/pow(inv_phMui,2));
        Rprintf("Samp %d inv_phMui^2: %f\n",i+1,pow(inv_phMui,2));*/
    }
    
    double score = score1 * (-inv_ph*inv_ph) + 2*lambda*ph;
    double info = info1 * pow(inv_ph,4) + 2*lambda;
    
    return List::create(score,info);
}

int sign(double x) {
    return (x > 0) - (x < 0);
}

double soft_thresh(double theta, double lambda, double alpha){ /* ST of SCAD penalty */
  double STval;
	double a=3.7;
	
    if(fabs(theta)<=2*lambda*alpha){
  		if(fabs(theta)<=lambda*alpha){
  			STval = 0;
  		} else{
  			STval = sign(theta)*(fabs(theta)-alpha/(1-alpha));
  		}
    } else if(fabs(theta)>2*lambda*alpha && fabs(theta)<=a*lambda*alpha){
        STval = ((a-1)*theta - sign(theta)*a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha)));
    } else{
		STval = theta;
	}
    return(STval);
}

double phi_ml(arma::vec y, arma::vec mu, arma::rowvec wts, int limit, int trace){
    double eps = 0.0001220703; /* from R */
    double p0 = 0;
    double N = accu(wts);
    int n = y.size();
    arma::vec wtd_y(n);
    
    /*arma::uvec wts_one_ids = find(wts == 1);
    arma::vec y_one = y(wts_one_ids);*/
    
    int equal=0;
    if(n==1){
        return(p0=0);
    } else {
        for(int i=0; i<n; i++){
            wtd_y(i)=y(i)*wts(i);
        }
    }
    
    double sum_wtd_y = accu(wtd_y);
    for(int i=0; i<n; i++){
        if(wtd_y(i)==sum_wtd_y){
            return(p0=0);
        }
        if(i<n-1){
            if(wtd_y(i) == wtd_y(i+1)){
                equal++;
            }
        }
    }
    
    if(equal==(n-1)){
        return(p0=0);
    }
    if(accu(mu)<1e-13){
        return(p0=0);
    }
    
    /*
    int all_unequal=0;
    if(y_one.size()==0){
        all_unequal=1;
    } else {
        for(int i=0; i<y_one.size(); i++){
            if(i > 0){
                if(y_one(i)!=y_one(i-1)){
                    all_unequal=1;
                }
            }
        }
    }
    
    if(all_unequal==0){
        return(p0=0);
    }
    */
    
    for(int i=0; i<n; i++){
        double wtsi=wts(i), yi=y(i), mui=mu(i);
        p0 += wtsi * pow(yi/mui-1,2);
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




/* */



// [[Rcpp::export]]
List M_step(int j, int a, arma::vec y_j, arma::mat all_wts, arma::vec offset, int k, arma::mat theta, arma::vec coefs_j, arma::vec phi_j, int cl_phi, double phi_g, int est_phi, double lambda, double alpha, double IRLS_tol, int maxit_IRLS){

    arma::vec beta = coefs_j;
    arma::mat temp(maxit_IRLS, (2*k));
    arma::mat temp_beta(maxit_IRLS,k);
    arma::mat temp_phi(maxit_IRLS,k);
    int continue_beta = 1;       /* Continue estimating beta as long as it's 1 (set to 0 when under SSE IRLS tol level)*/
    int continue_phi = 1;
    temp.zeros();
    
    
    int n = y_j.size();
    arma::mat eta(n,k), mu(n,k);
    eta.zeros();
    mu.zeros();
    
    /*
    arma::vec eta_g(n),mu_g(n);
    eta_g.zeros();
    mu_g.zeros();
    
    
    double logmeany = log(mean(y_j));
    
    for(int ii=0;ii<n;ii++){
        eta_g(ii) = logmeany + offset(ii);
        mu_g(ii) = exp(eta_g(ii));
        Rprintf("Sample %d: eta_g = %f, mu_g = %f",ii+1, eta_g(ii),mu_g(ii));
    } */
    
    /*arma::rowvec n_ones(n);
    n_ones.ones();
    
    if(cl_phi==0){
        phi_g = phi_ml(y_j,mu_g,n_ones,10,0);
        for(int c=0;c<k;c++){
            phi_j(c) = phi_g(j-1);
        }
    } */ /* THIS IS WRONG */
    
    
    /* IRLS */
    for(int i=0; i<maxit_IRLS; i++){
        /* Initiate eta */
        if(i==0) {
            for(int ii=0; ii<n; ii++) {
                for(int jj=0; jj<k; jj++) {
                    if(a==1){
                        eta(ii,jj) = beta(jj);
                    } else {
                        eta(ii,jj) = beta(jj) + offset(ii);
                    }
                    /*mu(ii,jj) = exp(eta(ii,jj));*/
                    mu(ii,jj) = pow(2,eta(ii,jj));
                }
            }
        }
        
        if(est_phi==1 && cl_phi==0 && continue_phi==1){
          phi_g = phi_ml_g(y_j,mu,all_wts,10,0);
          for(int c=0; c<k; c++){
            phi_j(c) = phi_g;
          }
        }
        
        
        
        /* Initiate temp matrix to track IRLS */
        int idx = 0;
        for(int ii=0; ii<k; ii++){
            temp(i,idx) = beta(ii);
            idx++;
        }
        for(int ii=0; ii<k; ii++){
            temp(i,idx) = phi_j(ii);
            idx++;
        }
        
        for(int ii=0; ii<k; ii++){
          temp_beta(i,ii) = beta(ii);
          temp_phi(i,ii) = phi_j(ii);
        }
    
        /* CDA */
        for(int c=0; c<k; c++){
            
            /* Testing gene-wise phi (ADD TO INPUT IN RCPP "arma::vec fixed_phi" AND FUNCTION & LogLike IN R AS WELL) */
            /* phi_j(c) = fixed_phi(j-1); */
    
            
            arma::rowvec wts_c = all_wts.row(c);
    
            /* First calculate all trans_y, w, and products */
            arma::vec all_trans_y(n);
            arma::vec all_w(n);
            arma::vec all_prod_w_trans_y(n);
    
            for(int ii=0; ii<n; ii++){
                all_trans_y(ii) = ( eta(ii,c)-offset(ii) ) + (y_j(ii)-mu(ii,c))/mu(ii, c);
                all_w(ii) = sqrt(wts_c(ii)*mu(ii,c)*mu(ii,c)/(mu(ii,c)+mu(ii,c)*mu(ii,c)*phi_j(c)));
                all_prod_w_trans_y(ii) = all_trans_y(ii)*all_w(ii);
            }
    
            /* Subset just the values of trans_y, w, and products where weight != 0 */
            arma::uvec good_ids = find(wts_c != 0);
    
            arma::vec trans_y = all_trans_y.rows(good_ids);
            arma::vec w = all_w.rows(good_ids);
            arma::vec prod_w_trans_y = all_prod_w_trans_y.rows(good_ids);
    
            /* Update beta */
            if(continue_beta==1){
              if((1-alpha)*lambda != 0){                    /* MLE update w/ SCAD*/
                  beta(c) = ((1-alpha)*lambda*((accu(beta)-beta(c))+(accu(theta.row(c))-theta(c,c))) + accu(prod_w_trans_y)/n )  / ((1-alpha)*lambda*(k-1) + accu(w)/n );
              } else {
                  beta(c) = accu(prod_w_trans_y)/accu(w);
              }
    
              if(beta(c) < (-30)){
                  /* Rprintf("Cluster %d, gene %d truncated at -30",c+1,j); */
                  beta(c) = -30;
              } else if(beta(c)>30){
                  /* Rprintf("Cluster %d, gene %d truncated at +30",c+1,j); */
                  beta(c) = 30;
              }
            }
    
            for(int ii=0; ii<n; ii++){
                eta(ii,c) = beta(c) + offset(ii);
                /*mu(ii,c) = exp(eta(ii,c));*/
                mu(ii,c) = pow(2,eta(ii,c));
            }
    
            /* Estimate phi */
            //Rprintf("phi.ml iter %d, cluster %d \n",i+1,c+1);
            if(cl_phi==1 && est_phi==1 && continue_phi==1){
              phi_j(c) = phi_ml(y_j,mu.col(c),wts_c,10,0);
            }
            
            /* Testing fixing phi to be = 10 */
            /* phi_j(c) = 10; */
            
        }
    
    
        /* Update theta matrix */
        for(int cc=0; cc<k; cc++){
            for(int ccc=0; ccc<k; ccc++){
                 theta(cc,ccc) = soft_thresh(beta(cc)-beta(ccc),lambda,alpha);
            }
        }
        
        if(i>0){
            double SSE_beta=0;
            double SSE_phi=0;
            for(int cc=0; cc<k; cc++){
              SSE_beta += pow(temp_beta(i,cc)-temp_beta(i-1,cc),2);
              SSE_phi += pow(temp_phi(i,cc)-temp_phi(i-1,cc),2);
            }
            Rprintf("SSE beta: %f, SSE phi: %f\n",SSE_beta,SSE_phi);
            if(SSE_beta<IRLS_tol){
              continue_beta=0;
            }
            if(SSE_phi<IRLS_tol){
              continue_phi=0;
            }
        }
        
        if(i==maxit_IRLS-1){
          break;
        }
        if(continue_beta==0 && continue_phi==0){
          break;
        }
    }
    
    
    return List::create(Rcpp::Named("coefs_j")=beta,Rcpp::Named("theta_j")=theta,Rcpp::Named("temp_j")=temp,Rcpp::Named("phi_j")=phi_j);
    
}
