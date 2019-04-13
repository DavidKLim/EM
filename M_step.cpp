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
#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

int sign(double x) {
    return (x > 0) - (x < 0);
}

double SCAD_soft_thresh(double theta, double lambda, double alpha){ /* ST of SCAD penalty */
  double STval;
	double a=3.7;
	
    /*if(fabs(theta)<=2*lambda*alpha){*/
	if(fabs(theta)<=(alpha/(1-alpha))+lambda*alpha){
  		if(fabs(theta)<=alpha/(1-alpha)){
  			STval = 0;
  		} else{
  			STval = sign(theta)*(fabs(theta)-alpha/(1-alpha));
  		}
    } else if(fabs(theta)>(alpha/(1-alpha))+lambda*alpha && fabs(theta)<=a*lambda*alpha){
		double omega = ((a-1)*theta)/(a-1-1/(lambda*(1-alpha)));
		if(fabs(omega) - (a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha))) <= 0){
			STval=0;
		} else{
			STval = sign(omega)*(fabs(omega)-(a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha))));
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

List score_info(int N, double ph, arma::vec mu, arma::vec y, arma::vec wts){
    double lambda = 1e-25;
    double score1 = 0, info1 = 0;
    double inv_ph = 1/ph + lambda;  // try adding lambda to stabilize (when phi very large)?
    ph=ph/(1+ph*1e-25);
    double scorei, infoi;
    
    int n = y.size();
	
	/*Timer timer2;
	timer2.step("Initial");*/
	
	/* Vectorizing and using Rcpp sugar digamma/trigamma functions: not really faster */
	/*NumericVector scoreic(n*k),infoic(n*k),yic(n*k),muic(n*k),wtsic(n*k),inv_phic(n*k);
	inv_phic.fill(inv_ph);
	int index=0;
	for(int i=0; i<n; i++){
        for(int c=0; c<k; c++){
			yic(index)=y(i);
			muic(index)=mu(i,c);
			wtsic(index)=wts(c,i);
			index++;
		}
	}
	timer2.step("Initialization");
	
	scoreic=wtsic * (digamma(inv_phic+yic) - digamma(inv_phic) + log(inv_phic) + 1 - log(inv_phic+muic)-(yic+inv_phic)/(inv_phic+muic));
	infoic=wtsic * (trigamma(inv_phic) + 2/(inv_phic+muic) - trigamma(inv_phic+yic) - ph - (yic+inv_phic)/pow(inv_phic+muic,2));
	
	timer2.step("Score/info calculation");
	score1=sum(scoreic);
	info1=sum(infoic);
	timer2.step("Summation");
	NumericVector res2(timer2);
	Rcpp::print(res2);*/
	
	for(int i=0; i<n; i++){			
            
        scorei = wts(i) * (R::digamma(inv_ph+y(i)) - R::digamma(inv_ph) + log(inv_ph) + 1 - log(inv_ph+mu(i)) - (y(i)+inv_ph)/(inv_ph+mu(i)));
        infoi = wts(i) * (R::trigamma(inv_ph) + 2/(inv_ph+mu(i)) - R::trigamma(inv_ph+y(i)) - ph - (y(i)+inv_ph)/pow(inv_ph+mu(i),2));
            
        score1 += scorei;
        info1 += infoi;
	}
	
	/*timer2.step("scoreic/infoic calc");
	NumericVector res2(timer2);
	Rcpp::print(res2);*/
	
    
    double score = score1 * (-inv_ph*inv_ph) + 2*lambda*ph;
    double info = info1 * pow(inv_ph,4) + 2*lambda;
    
    return List::create(score,info);
}

double phi_ml_g(arma::vec y, arma::vec mu, arma::vec wts, int limit, int trace){
    double eps = 0.0001220703; /* from R */
    double p0 = 0;
    double N = accu(wts);
    int n = y.size();
    
    if(n==1){
        return(p0=0);
    }
    if(accu(mu)<1e-13){
        return(p0=0);
    }
    
	for(int i=0; i<n; i++){
		p0 += wts(i)*pow(y(i)/mu(i)-1,2);
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
	if(trace==1){Rprintf("phi_ml() converged after %d iterations\n",it);}
    
    return(p0);
}

double phi_ml(arma::vec y, arma::vec mu, arma::vec wts, int limit, int trace){
    double eps = 0.0001220703; /* from R */
    double p0 = 0;
    double N = accu(wts);
    int n = y.size();
    arma::vec wtd_y(n);
    
    
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
            return(p0=0);		/* Case where there is only one non-zero sample in a cluster */
        }
        if(i<n-1){
            if(wtd_y(i) == wtd_y(i+1)){
                equal++;        /* Tracking how many samples have equal expression */
            }
        }
    }
    
    if(equal==(n-1)){
        return(p0=0);           /* If all samples are equal in expression --> phi =0 */
    }
    if(accu(mu)<1e-13){
        return(p0=0); 
    }
    
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
List M_step(arma::mat X, int p, int j, int a, arma::vec y_j, arma::mat all_wts, arma::vec vec_wts, arma::vec keep, arma::vec offset, int k, arma::mat theta, arma::vec coefs_j, arma::vec phi_j, int cl_phi, double phi_g, int est_phi, int est_covar, double lambda, double alpha, double IRLS_tol, int maxit_IRLS){

    arma::vec beta = coefs_j;
    arma::mat temp(maxit_IRLS, (2*k+p));
    arma::mat temp_beta(maxit_IRLS,k+p);
    arma::mat temp_phi(maxit_IRLS,k);
    int continue_beta = 1;       /* Continue estimating beta as long as it's 1 (set to 0 when under SSE IRLS tol level)*/
    int continue_phi = 1;
	int continue_gamma = 1;
    temp.zeros();
    
    int n = y_j.size();
	arma::vec n_k(k);
	arma::vec rep_y_j(n*k);
	for(int c=0; c<k; c++){
		n_k(c) = sum(all_wts.row(c));
		for(int i=0; i<n; i++){
			rep_y_j(c*n+i) = y_j(i);
		}
	}
	
	int index=0;
	
	arma::vec y_tilde(n*k);
	arma::vec eta(n*k), mu(n*k);
	arma::mat mat_mu(n,k);
	/*arma::mat mat_W(n*k,n*k);
	mat_W.zeros();*/
	arma::vec vec_W(n*k);
	vec_W.zeros();
	arma::vec resid(n*k);
	
	/* Turn on this for WLS estimate of covariates */
	arma::vec MLE_beta(p+k);
	arma::vec beta_cls(k);
	
	Timer timer;
	int output_timer = 0;
	int TIME_ITER = 0;
	
	/* PP filtering ids, gene */
	arma::uvec ids = find(keep == 1);
	/*arma::uvec ids = find(keep >= 0);*/   /* test: keep all samples */
	
	arma::uvec ids_c = find(keep==5);     /* just dummy initialization */

	
    /* IRLS */
    for(int i=0; i<maxit_IRLS; i++){
		if(i==TIME_ITER){
			timer.step("start");
		}
		
		/* Calculate eta and mu */
		eta = X * beta + offset;
		
		if(i==TIME_ITER){
			timer.step("eta");
		}
		
		for(int ii=0; ii<(n*k); ii++){
			mu(ii)=pow(2,eta(ii));
			mat_mu(ii-n*floor(ii/n),floor(ii/n)) = mu(ii);
		}
		
		if(i==TIME_ITER){
			timer.step("mu/mat_mu");
		}

        
        /* Initiate temp matrix to track IRLS */
        int idx = 0;
        for(int c=0; c<k+p; c++){
            temp(i,idx) = beta(c);
			temp_beta(i,c) = beta(c);
            idx++;
        }
        for(int c=0; c<k; c++){
            temp(i,idx) = phi_j(c);
			temp_phi(i,c) = phi_j(c);
			beta_cls(c) = beta(c);
            idx++;
        }
		
		if(i==TIME_ITER){
			timer.step("temp mats");
		}
		
		/* Estimating phi */
		if(phi_g == 0 || a>1 || i>0){              /* don't estimate phi if glm.nb() was fit to initialize phi */
			if(est_phi==1 && cl_phi==0 && continue_phi==1){
			  phi_g = phi_ml_g(rep_y_j(ids),mu(ids),vec_wts(ids),10,0);
			  if(phi_g>100){
					phi_g=100;
				}
			  for(int c=0; c<k; c++){
				phi_j(c) = phi_g;
			  }
			}
		}
		/*Rprintf("phi iter%d: %f\n",i,phi_g);*/
		
		if(i==TIME_ITER){
			timer.step("est phi_g");
		}
		
		
		/* Calculate y_tilde and W matrix */
		for(int cc=0; cc<k; cc++){
			for(int ii=0; ii<n; ii++){
				index=ii+(n*cc);
				y_tilde(index) = ( eta(index)-offset(ii) ) + (y_j(ii)-mu(index))/(mu(index)*log(2));      /* link function: log_2(mu) = eta */
				double w_ii = sqrt(all_wts(cc,ii)*mu(index)*mu(index)/(mu(index)+mu(index)*mu(index)*phi_j(cc)));
				/*mat_W(index,index)=w_ii;*/    /* mat_W need not be calculated */
				vec_W(index)=w_ii;
			}
		}
		
		if(i==TIME_ITER){
			timer.step("calc y_tilde/vec_W");
		}
		
		/* Hard coded coordinate-wise covar_beta (something is wrong here): */
		if(est_covar==1 && continue_gamma==1){
			
			/* WLS. Alternate: coordinate-wise covariate updates */
			/*MLE_beta = inv(X.t() * mat_W * X) * X.t() * mat_W * y_tilde;
			for(int pp=k; pp<(p+k); pp++){
				beta(pp)=MLE_beta(pp);
			}*/
			
			/* Hard coded coordinate-wise covar_beta: */
			for(int pp=k; pp<(p+k); pp++){
				arma::mat Xpp=X;
				arma::vec betapp=beta;
				Xpp.shed_col(pp);
				betapp.shed_row(pp);			
				resid=y_tilde-Xpp*betapp;
				
				arma::vec subs_vec_W=vec_W(ids);
				arma::vec Xcolpp=X.col(pp);
				arma::vec subs_Xcolpp=Xcolpp(ids);
				arma::vec subs_resid=resid(ids);             /* subsetting just PP >0.001 */
				
				beta(pp) = accu(subs_vec_W % subs_Xcolpp % subs_resid)/accu(subs_vec_W % pow(subs_Xcolpp,2));
				
				/*Rprintf("covar iter%d gamma=%f, top=%f, bottom=%f\n",i,beta(pp),accu(subs_vec_W % subs_Xcolpp % subs_resid),accu(subs_vec_W % pow(subs_Xcolpp,2)));*/
				if(beta(pp)>100){
					beta(pp)=100;
				} else if(beta(pp)<-100){
					beta(pp)=-100;
				}   /* is this feasible to do? 2^100 is 1E30, 2^50 is 1.1E15 */
			}
		}
		
		if(i==TIME_ITER){
			timer.step("cov ests");
		}
		
        /* CDA */
        for(int c=0; c<k; c++){
			beta_cls(c) = beta(c);
			arma::uvec fused_ids = find(theta.row(c) == 0);        /* fused cluster labels */
			arma::uvec not_fused_ids = find(theta.row(c) != 0);
			int fused_n_k=accu(n_k(fused_ids));                   /* n_k = number of samples in fused cluster (if fused) */
			int num_fused_cls = fused_ids.n_elem;                 /* tracks number of clusters that are fused with current cl c */
			/*Rprintf("%d\n",fused_ids.n_elem);*/

			if(num_fused_cls<=1){
				ids_c = find(X.col(c) % keep == 1);             /* if no fused cluster, subset to just PP > 0.001 samples in cl c */
				/*if(c==0){ids_c.print("ids2");}*/
			} else{
				int min_fused_id = fused_ids.index_min();       /* determine which of the fused clusters has the smallest label */
				int min_fused_cl = fused_ids(min_fused_id);     
				/*Rprintf("min:%d,cl%d",min_fused_cl,c);*/
				if(c > min_fused_cl){                           /* if current cl c is fused w/ a previous cl, then take the prev ests. just beta? */
					beta(c) = beta(min_fused_cl); 
					/*if(cl_phi==1 && est_phi==1 && continue_phi==1){
						phi_j(c) = phi_j(min_fused_cl);
					}*/
					/*Rprintf("Set beta%d and phi_j%d equal to beta%d and phi_j%d\n",c,c,min_fused_cl,min_fused_cl);*/
				}
				arma::mat arma_Xsum = X.cols(fused_ids);
				NumericMatrix Xsum = wrap(arma_Xsum);
				NumericVector Xfused = rowSums(Xsum);
				arma::vec arma_Xfused = as<arma::vec>(Xfused);      /* sum manipulations: essentially combine X matrix of fused cls by rowSumming */
				ids_c = find(arma_Xfused % keep == 1);             /* PP filt >0.001 on new sample size of combined cls */
				/*if(c==0){ids_c.print("ids1");}*/
			}
			arma::vec resid_c(ids_c.n_elem);
            
			/*if(c==0){ids_c.print();}*/
			
			if(p==0){
				resid_c=y_tilde(ids_c);       		/* if no covars: resid = y-Xb^(-1) = y */
			} else{
				arma::mat Xc=X.rows(ids_c);
				arma::vec betac=beta;
				if(num_fused_cls<=1){
					Xc.shed_col(c);
					betac.shed_row(c);
				} else{                             /* if fused clusters, then remove in Xc*betac both cluster samples */
					Xc=Xc.cols(not_fused_ids);
					betac=betac.rows(not_fused_ids);
				}
				resid_c=y_tilde(ids_c)-Xc*betac;	
			}

			
			arma::vec vec_W_c = vec_W(ids_c);
			
            /* Update beta */
            if(continue_beta==1){
				
			  if((1-alpha)*lambda != 0){
                  beta(c) = ((1-alpha)*lambda*((accu(beta_cls)-beta(c))+(accu(theta.row(c))-theta(c,c))) + accu(vec_W_c % resid_c)/(fused_n_k) )  /
				  ((1-alpha)*lambda*(k-1) + accu(vec_W_c)/(fused_n_k) );
              } else {
				  beta(c) = accu(vec_W_c % resid_c)/accu(vec_W_c);
              }
			  
				  /*Rprintf("Iter%d cl%d:%f, 1=%f, 2=%f, 3=%f, 4=%f\n",i,c,beta(c),
				  (1-alpha)*lambda*((accu(beta_cls)-beta(c))+(accu(theta.row(c))-theta(c,c))),accu(vec_W_c % resid_c)/(n_k(c)),(1-alpha)*lambda*(k-1),accu(vec_W_c)/(n_k(c)));*/
    
              if(beta(c) < (-100)){
                  /* Rprintf("Cluster %d, gene %d truncated at -100",c+1,j); */
                  beta(c) = -100;
              } else if(beta(c)>100){
                  /* Rprintf("Cluster %d, gene %d truncated at +100",c+1,j); */
                  beta(c) = 100;
              }
            }

			/*eta = X * beta + offset;
            for(int ii=0; ii<(n*k); ii++){
                mu(ii) = pow(2,eta(ii));
            }*/
    
            /* Estimate phi */
            /*Rprintf("phi.ml iter %d, cluster %d \n",i,c);*/
            if(cl_phi==1 && est_phi==1 && continue_phi==1){
				phi_j(c) = phi_ml(rep_y_j(ids_c),mu(ids_c),vec_wts(ids_c),10,0);
            }
			
            
        }
		
		if(i==TIME_ITER){
			timer.step("CDA on log2 baselines");
		}
		
			
		/* Update theta matrix */
		for(int cc=0; cc<k; cc++){
			for(int ccc=0; ccc<k; ccc++){
				theta(cc,ccc) = SCAD_soft_thresh(beta(cc)-beta(ccc),lambda,alpha);
				/*theta(cc,ccc) = lasso_soft_thresh(beta(cc)-beta(ccc),lambda*alpha);*/
			}
		}

		if(i==TIME_ITER){
			timer.step("Theta");
		}
        
        if(i>2){
            double diff_beta=0;
			double diff_gamma=0;
            double diff_phi=0;
            for(int cc=0; cc<k; cc++){
              diff_beta += fabs(temp_beta(i,cc)-temp_beta(i-1,cc))/fabs(temp_beta(i-1,cc)*k);
			  if(cl_phi==1){
				diff_phi += fabs(temp_phi(i,cc)-temp_phi(i-1,cc))/fabs(temp_phi(i-1,cc)*k);
			  } else{
				diff_phi = fabs(temp_phi(i,cc)-temp_phi(i-1,cc))/fabs(temp_phi(i-1,cc));
			  }
            }
			for(int cc=k; cc<(k+p); cc++){
				if(p>0){
					diff_gamma += fabs(temp_beta(i,cc)-temp_beta(i-1,cc))/fabs(temp_beta(i-1,cc)*p);
				}
			}
            /*Rprintf("diff beta: %f, diff phi: %f\n, diff gamma: %f\n",diff_beta,diff_phi,diff_gamma);*/
			
            if(diff_beta<IRLS_tol){
              continue_beta=0;
            }
			if(p>0){
				if(diff_gamma<IRLS_tol){
					continue_gamma=0;
				}
			} else{continue_gamma=0;}
            if(diff_phi<IRLS_tol){
              continue_phi=0;
            }
        }
        
        if(i==maxit_IRLS-1){
          break;
        }
        if(continue_beta==0 && continue_gamma==0 && continue_phi==0){
          break;
        }
		
		if(i==TIME_ITER){
			timer.step("break conds");
		}
		
		if(i==TIME_ITER && output_timer==1){
			NumericVector res(timer);
			Rcpp::print(res);
		}
    }
    
    
    return List::create(Rcpp::Named("coefs_j")=beta,Rcpp::Named("theta_j")=theta,Rcpp::Named("temp_j")=temp,Rcpp::Named("phi_j")=phi_j);
    
}
