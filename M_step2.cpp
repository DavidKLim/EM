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

double SCAD_soft_thresh(double theta, double lambda, double alpha){ /* ST of SCAD penalty */
  double STval;
	double a=3.7;
	
    if(fabs(theta)<=2*lambda*alpha){
  		if(fabs(theta)<=lambda*alpha){
  			STval = 0;
  		} else{
  			STval = sign(theta)*(fabs(theta)-alpha/(1-alpha));
  		}
    } else if(fabs(theta)>2*lambda*alpha && fabs(theta)<=a*lambda*alpha){
		double mid_term = ((a-1)*theta - a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha)));
		if(mid_term > 0){
			STval = sign( ((a-1)*theta)/(a-1-1/(lambda*(1-alpha))) )*mid_term;
		} else{
			STval=0;
		}
    } else{
		STval = theta;
	}
    return(STval);
}

double lasso_soft_thresh(double alpha, double lambda){
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
List M_step(arma::mat X, int p, int j, int a, arma::vec y_j, arma::mat all_wts, arma::vec offset, int k, arma::mat theta, arma::vec coefs_j, arma::vec phi_j, int cl_phi, double phi_g, int est_phi, int est_covar, double lambda, double alpha, double IRLS_tol, int maxit_IRLS){

    arma::vec beta = coefs_j;
    arma::mat temp(maxit_IRLS, (2*k+p));
    arma::mat temp_beta(maxit_IRLS,k+p);
    arma::mat temp_phi(maxit_IRLS,k);
    int continue_beta = 1;       /* Continue estimating beta as long as it's 1 (set to 0 when under SSE IRLS tol level)*/
    int continue_phi = 1;
    temp.zeros();
    
    
    int n = y_j.size();
	arma::vec n_k(k);
	for(int c=0; c<k; c++){
		n_k(c) = sum(all_wts.row(c));
	}
	
	int index=0;
	
	arma::vec y_tilde(n*k);
	arma::vec eta(n*k), mu(n*k);
	arma::mat mat_mu(n,k);
	arma::mat mat_W(n*k,n*k);
	mat_W.zeros();
	arma::vec vec_W(n*k);
	vec_W.zeros();
	arma::vec resid(n*k);
	
	arma::vec vec_wts(n*k);
	for(int i=0; i<n; i++){
		for(int c=0; c<k; c++){
			index=i+(n*c);
			vec_wts(index) = all_wts(c,i);
		}
	}
	
	/* Turn on this for WLS estimate of covariates */
	arma::vec MLE_beta(p+k);
	arma::vec beta_cls(k);
	
	
    /* IRLS */
    for(int i=0; i<maxit_IRLS; i++){
		
		/* Calculate eta and mu */
		eta = X * beta + offset;
		for(int ii=0; ii<(n*k); ii++){
			mu(ii)=pow(2,eta(ii));
		}
		
		/* Calculate working response y_tilde and mat_mu */
		for(int c=0; c<k; c++){
			for(int ii=0; ii<n; ii++){
				index=ii+(n*c);
				y_tilde(index) = ( eta(index)-offset(ii) ) + (y_j(ii)-mu(index))/(mu(index)*log(2));      /* link function: log_2(mu) = eta */
				mat_mu(ii,c) = mu(index);                                                                 /* d(log2(mu))/d(mu) = 1/(mu*ln(2)) */
			}
		}		
		
		/* Estimating phi */
        if(est_phi==1 && cl_phi==0 && continue_phi==1){
          phi_g = phi_ml_g(y_j,mat_mu,all_wts,10,0);
		  /*if(phi_g>5){
			  phi_g=5;
		  }*/
          for(int c=0; c<k; c++){
            phi_j(c) = phi_g;
          }
        }
		
		/* Calculate W matrix */
		for(int cc=0; cc<k; cc++){
			for(int ii=0; ii<n; ii++){
				index=ii+(n*cc);
				mat_W(index,index)=sqrt(all_wts(cc,ii)*mu(index)*mu(index)/(mu(index)+mu(index)*mu(index)*phi_j(cc)));
				vec_W(index)=mat_W(index,index);
			}
		}
        
		/*Rprintf("IRLS iter %d\n phi: %f\n",i,phi_g);*/
		
		/* WLS. Alternate: hard-coding covariate updates */
		if(est_covar==1){
			MLE_beta = inv(X.t() * mat_W * X) * X.t() * mat_W * y_tilde;
			for(int pp=k; pp<(p+k); pp++){
				beta(pp)=MLE_beta(pp);
			}
		}
		
		/* */
		/* Hard coded covar_beta (something is wrong here): */
			/*for(int pp=k; pp<(p+k); pp++){
				arma::mat Xpp=X;
				arma::vec betapp=beta;
				Xpp.shed_col(pp);
				betapp.shed_row(pp);			
				resid=y_tilde-Xpp*betapp;
				
				beta(pp) = accu(vec_W % X.col(p) % resid)/accu(vec_W % pow(X.col(p),2));
		}*/
        
        /* Initiate temp matrix to track IRLS */
        int idx = 0;
        for(int c=0; c<k+p; c++){
            temp(i,idx) = beta(c);
            idx++;
        }
        for(int c=0; c<k; c++){
            temp(i,idx) = phi_j(c);
            idx++;
        }
		for(int c=0; c<k+p; c++){
			temp_beta(i,c) = beta(c);
		}
        for(int c=0; c<k; c++){
          temp_phi(i,c) = phi_j(c);
        }
		for(int c=0; c<k; c++){
			beta_cls(c) = beta(c);
		}
		
        /* CDA */
        for(int c=0; c<k; c++){
            
            /* Testing gene-wise phi (ADD TO INPUT IN RCPP "arma::vec fixed_phi" AND FUNCTION & LogLike IN R AS WELL) */
            /* phi_j(c) = fixed_phi(j-1); */
			arma::mat Xc=X;
			arma::vec betac=beta;
			Xc.shed_col(c);
			betac.shed_row(c);
			resid=y_tilde-Xc*betac;

			beta_cls(c) = beta(c);
			
			
			
            /* Update beta */
            if(continue_beta==1){
              if((1-alpha)*lambda != 0){                    /* MLE update w/ SCAD*/
                  beta(c) = ((1-alpha)*lambda*((accu(beta_cls)-beta(c))+(accu(theta.row(c))-theta(c,c))) + accu(vec_W % X.col(c) % resid)/(n*k) )  /
				  ((1-alpha)*lambda*(k-1) + accu(vec_W % pow(X.col(c),2)/(n*k) ));
              } else {
				  beta(c) = accu(vec_W % X.col(c) % resid)/accu(vec_W % pow(X.col(c),2));
				  /*Rprintf("cl %d::\n",c);
				  Rprintf("top: %f\n",accu(vec_W % X.col(c) % resid));
				  Rprintf("bottom: %f\n",accu(vec_W % pow(X.col(c),2)));
				  vec_W.t().print("weights:");
				  X.col(c).t().print("X.col(c):");
				  Rprintf("beta: %f\n",beta(c));*/
              }
    
              if(beta(c) < (-30)){
                  /* Rprintf("Cluster %d, gene %d truncated at -30",c+1,j); */
                  beta(c) = -30;
              } else if(beta(c)>30){
                  /* Rprintf("Cluster %d, gene %d truncated at +30",c+1,j); */
                  beta(c) = 30;
              }
            }

			eta = X * beta + offset;
            for(int ii=0; ii<(n*k); ii++){
                mu(ii) = pow(2,eta(ii));
            }
    
            /* Estimate phi */
            //Rprintf("phi.ml iter %d, cluster %d \n",i+1,c+1);
            if(cl_phi==1 && est_phi==1 && continue_phi==1){
              phi_j(c) = phi_ml(y_j,mat_mu.col(c),all_wts.row(c),10,0);
            }
			
			/* Update theta matrix */
			for(int cc=0; cc<k; cc++){
				for(int ccc=0; ccc<k; ccc++){
					theta(cc,ccc) = SCAD_soft_thresh(beta(cc)-beta(ccc),lambda,alpha);
					/*theta(cc,ccc) = lasso_soft_thresh(beta(cc)-beta(ccc),lambda*alpha);*/
				}
			}
            
            /* Testing fixing phi to be = 10 */
            /* phi_j(c) = 10; */
            
        }
		
		/* Calculate working response y_tilde and mat_mu */
		/*for(int c=0; c<k; c++){
			for(int ii=0; ii<n; ii++){
				index=ii+(n*c);
				y_tilde(index) = ( eta(index)-offset(ii) ) + (y_j(ii)-mu(index))/mu(index);
				mat_mu(ii,c) = mu(index);
			}
		}		*/
		/* Calculate W matrix */
		/*for(int cc=0; cc<k; cc++){
			for(int ii=0; ii<n; ii++){
				index=ii+(n*cc);
				mat_W(index,index)=sqrt(all_wts(cc,ii)*mu(index)*mu(index)/(mu(index)+mu(index)*mu(index)*phi_j(cc)));
				vec_W(index)=mat_W(index,index);
			}
		}*/
		
			



        
        if(i>1){
            double SSE_beta=0;
            double SSE_phi=0;
            for(int cc=0; cc<k; cc++){
              SSE_beta += pow(temp_beta(i,cc)-temp_beta(i-1,cc),2);
              SSE_phi += pow(temp_phi(i,cc)-temp_phi(i-1,cc),2);
            }
            /* Rprintf("SSE beta: %f, SSE phi: %f\n",SSE_beta,SSE_phi); */
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
