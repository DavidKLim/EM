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

/* Create functions phi_ml() and soft_thresholding() here */

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

double soft_thresh(double alpha, double lambda){
    double STval;
    if(fabs(alpha)-lambda<0){
        STval = 0;
    } else {
        STval = sign(alpha) * (fabs(alpha)-lambda);
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
List M_step(int j, int a, arma::vec y_j, arma::mat all_wts, arma::vec offset, int k, List theta_list, arma::mat coefs, arma::mat phi, int cl_phi, arma::vec phi_g, double lambda1, double lambda2, double tau, double IRLS_tol, int maxit_IRLS){

    arma::mat beta = coefs.row(j-1), theta = theta_list[j-1], temp(maxit_IRLS, (2*k));
    temp.zeros();
    
    
    int n = y_j.size();
    arma::mat eta(n,k), mu(n,k);
    eta.zeros();
    mu.zeros();
    
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
                    mu(ii,jj) = exp(eta(ii,jj));
                }
            }
        }
        /* Initiate temp matrix to track IRLS */
        temp.row(i) = join_rows(beta.row(0),phi.row(j-1));
    
        /* CDA */
        for(int c=0; c<k; c++){
            
            /* Testing gene-wise phi (ADD TO INPUT IN RCPP "arma::vec fixed_phi" AND FUNCTION & LogLike IN R AS WELL) */
            /* phi(j-1,c) = fixed_phi(j-1); */
    
            
            arma::rowvec wts_c = all_wts.row(c);
    
            /* First calculate all trans_y, w, and products */
            arma::vec all_trans_y(n);
            arma::vec all_w(n);
            arma::vec all_prod_w_trans_y(n);
    
            for(int ii=0; ii<n; ii++){
                all_trans_y(ii) = ( eta(ii,c)-offset(ii) ) + (y_j(ii)-mu(ii,c))/mu(ii, c);
                all_w(ii) = sqrt(wts_c(ii)*mu(ii,c)*mu(ii,c)/(mu(ii,c)+mu(ii,c)*mu(ii,c)*phi(j-1,c)));
                all_prod_w_trans_y(ii) = all_trans_y(ii)*all_w(ii);
            }
    
            /* Subset just the values of trans_y, w, and products where weight != 0 */
            arma::uvec good_ids = find(wts_c != 0);
    
            arma::vec trans_y = all_trans_y.rows(good_ids);
            arma::vec w = all_w.rows(good_ids);
            arma::vec prod_w_trans_y = all_prod_w_trans_y.rows(good_ids);
    
            /* Update beta */
            if(lambda1 != 0){
                beta(c) = (lambda1*((accu(beta)-beta(c))+(accu(theta.row(c))-theta(c,c))) + accu(prod_w_trans_y)/n )  / (lambda1*(k-1) + accu(w)/n );
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
    
            for(int ii=0; ii<n; ii++){
                eta(ii,c) = beta(c) + offset(ii);
                mu(ii,c) = exp(eta(ii,c));
            }
    
            /* Estimate phi */
            //Rprintf("phi.ml iter %d, cluster %d \n",i+1,c+1);
            if(cl_phi==1){
                phi(j-1,c) = phi_ml(y_j,mu.col(c),wts_c,10,0);
            } else if(cl_phi==0){
                phi(j-1,c) = phi_g(j-1);
            }
            
            /* Testing fixing phi to be = 10 */
            /* phi(j-1,c) = 10; */
            
        }
    
    
        /* Update theta matrix */
        for(int cc=0; cc<k; cc++){
            for(int ccc=0; ccc<k; ccc++){
                if(fabs(theta(cc,ccc))>=tau){
                    theta(cc,ccc) = beta(cc)-beta(ccc);
                } else {
                    theta(cc,ccc) = soft_thresh(beta(cc)-beta(ccc),lambda2);
                }
            }
        }
        
        if(i>0){
            double SSE=0;
            for(int cc=0; cc<2*k; cc++){
                SSE += pow(temp(i,cc)-temp(i-1,cc),2);
            }
            if(SSE<IRLS_tol){
                break;
            }
        }
        if(i==maxit_IRLS-1){
            break;
        }
    }
    
    arma::mat final_phi = phi.row(j-1);
    
    return List::create(Rcpp::Named("coefs_j")=beta,Rcpp::Named("theta_j")=theta,Rcpp::Named("temp_j")=temp,Rcpp::Named("phi_j")=final_phi);
    
}
