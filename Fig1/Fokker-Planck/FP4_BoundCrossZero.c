
/*
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  File: FP4.c
;;;  Author: Roozbeh Kiani
;;;  Description: Fokker-Planck numerical solution, for temporal
;;;               propagation of a probability density based on a
;;;               diffusion process. allows for variable bound heights
;;;  Creation Date: Aug 15, 2006
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
*/



#include<stdio.h>
#include<math.h>

#include"mex.h"



#define LINEAREQ_TD_AllocMemX			0x10

long linearEq_Tridiagonal ( double *A , unsigned long EqNum , double *B , long flags , double **X )
 {
  double *backsubBuf;
  double *dbpA , *dbpB , *dbpX , *dbpBuf , divider;
  unsigned long i;

  if ( !A || !EqNum || !B )		return 0;
  if ( (backsubBuf=(double *)mxCalloc(EqNum,sizeof(double))) == NULL )		return 0;
  if ( flags & LINEAREQ_TD_AllocMemX )
   if ( (*X=(double *)mxCalloc(EqNum,sizeof(double))) == NULL )				return 0;

  dbpA = A;
  dbpB = B;
  dbpX = *X;
  divider=*(dbpA+1);
  *dbpX = (*dbpB)/divider;
  for ( i=1,dbpBuf=backsubBuf+1 ; i < EqNum ; i++,dbpBuf++ )
   { *dbpBuf = *(dbpA+2)/divider;
     dbpA += 3;
     divider = (*(dbpA+1)) - (*dbpA)*(*dbpBuf);
     dbpX++;
     dbpB++;
     *dbpX = ((*dbpB)-(*dbpA)*(*(dbpX-1)))/divider;		}
  for ( i=EqNum-1,dbpX=(*X)+EqNum-1,dbpBuf=backsubBuf+EqNum-1 ; i > 0 ; i--,dbpX--,dbpBuf-- )
   *(dbpX-1) -= (*dbpX)*(*dbpBuf);

  mxFree ( backsubBuf );
  return 1;
 }


/*
solving the Fokker-Planck equation for the diffusion process

Chang-Cooper method, a fully-implicit method
	see Chang and Cooper 1970, also Park and Petrosian

inputs:
	center_x		center of the dicrete spatial grid (corresponding to zero)
   dx				resolution of the dicrete spatial grid
   lb_margin	initial lower bound height
   ub_margin	initial upper bound height
   lb_change	change in the lower bound height in time, the length of this vector should be equal to stepnum
   ub_change	change in the upper bound height in time, the length of this vector should be equal to stepnum
   stepnum		number of steps (in time) that the probability distribution should be propagated
   dt				resolution of the discrete steps in time
   v           drift rate of the diffusion process
   sigma			standard deviation of the browninan motion noise

outputs:
   ufinal
   Pt				is the survivor function. it is the probability of NOT passing the bounds upto any moment
   Ptlow			probability of passing the lower bound at any moment
   Pthigh      probability of passing the upper bound at any moment
   Pg0			probability of being greater than zero at any moment
   Pxt			probability density across space and time
   
*/
long pde4 ( double *uinit , unsigned long ulen , unsigned long center_x , double dx ,
            double lb_margin , double ub_margin , double *lb_change , double *ub_change ,
            unsigned long stepnum , double dt , double v , double sigma ,
            double *ufinal , double *Pt , double *Ptlow , double *Pthigh , double *Pg0 , double *Pxt )
 {
  unsigned long x;
  double *ubuf;
  double *A;
  double Am , Bm , Cm , w , W , *Wplus , *Wminus;
  double *dbp2 , *dbp3;
  unsigned long t , r;
  long margin;

  if ( uinit==0 || ulen<2 || stepnum==0 || dx<=0 || dt<=0 )			return 0;

  		/*quantize lb_margin, ub_margin, lb_change and ub_change so that they refer to the indices on the
        xmesh rather than the actual values*/
  lb_margin = floor(lb_margin/dx+0.5);
  ub_margin = floor(ub_margin/dx+0.5);
  for ( t=0,dbp2=lb_change,dbp3=ub_change ; t < stepnum ; t++,dbp2++,dbp3++ )
   { *dbp2 = floor((*dbp2)/dx+0.5);
     *dbp3 = floor((*dbp3)/dx+0.5);	}

/*
  mexPrintf("ulen: %d\n center_x: %d\n dx: %f\n lb_margin: %f\n ub_margin: %f\n stepnum: %d\n dt: %f\n v: %f\n sigma: %f\n",
            ulen,center_x,dx,lb_margin,ub_margin,stepnum,dt,v,sigma);
*/

  		/*initialize the parameters of the algorithm*/
  Am = 1.0;
  Cm = 0.5*sigma*sigma;

  		/*allocate memory for the arrays used by the algorithm*/
  ubuf = (double *)mxCalloc(ulen+10,sizeof(double));
  Wplus = (double *)mxCalloc(ulen+10,sizeof(double));
  Wminus = (double *)mxCalloc(ulen+10,sizeof(double));
  A = (double *)mxCalloc(3*ulen,sizeof(double));

  		/*find the initial probability of having a value above zero*/
  for ( x=center_x+1,dbp2=Pg0,*dbp2=0,dbp3=uinit+center_x+1 ; x < ulen ; x++,dbp3++ )
	*dbp2 += fabs(*dbp3);
  *dbp2 += fabs(uinit[center_x])/2;

      /*find the probability of crossing the lower bound at time zero*/
  if ( lb_margin )
   for ( x=0,dbp2=Ptlow,*dbp2=0,dbp3=uinit; x < lb_margin ; x++,dbp3++ )
    { *dbp2 += fabs(*dbp3);		*dbp3 = 0;	}

    	/*find the probability of crossing the upper bound at time zero*/
  if ( ub_margin )
   for ( x=0,dbp2=Pthigh,*dbp2=0,dbp3=uinit+ulen-(long)ub_margin ; x < ub_margin ; x++,dbp3++ )
    { *dbp2 += fabs(*dbp3);		*dbp3 = 0;	}

  		/*find the initial probability of having a value between the two bounds*/
  for ( x=0,dbp2=Pt,*dbp2=0,dbp3=uinit ; x < ulen ; x++,dbp3++ )
  	*dbp2 += fabs(*dbp3);

  memcpy ( ufinal , uinit , ulen*sizeof(double) );

  		/*initialize the tri-diagonal matrix used by the Chang-Cooper fully-implicit method*/
  for ( x = 0 ; x <= ulen ; x++ )
   { Bm = -v;
     w = Bm*dx/Cm;
     if ( w == 0 )	W = 1.0;
     else				W = (w/2)/sinh(w/2);
     Wplus[x] = W * exp(w/2);
     Wminus[x] = W * exp(-w/2);	}
  for ( x = 0 ; x < ulen ; x++ )
   { A[3*x] = dt*Cm/(Am*dx*dx) * Wminus[x];
     A[3*x+1] = 1.0 + dt/(Am*dx)*( Cm/dx*Wplus[x] + Cm/dx*Wminus[x+1] );
     A[3*x+2] = dt*Cm/(Am*dx*dx) * Wplus[x+1];	}
  A[0] = 0.0;
  A[3*ulen-1] = 0.0;

  		/*now propagate the probability in time*/
  for ( t=0,r=1 ; t < stepnum ; t++,r++ )
   { memcpy ( ubuf , ufinal , ulen*sizeof(double) );
     linearEq_Tridiagonal ( A , ulen , ubuf , 0 , &ufinal );

     	/*find the probability of having a value greater than zero at this moment of time*/
     for ( x=center_x+1,dbp2=Pg0+r,*dbp2=0,dbp3=ufinal+center_x+1 ; x < ulen ; x++,dbp3++ )
     	*dbp2 += fabs(*dbp3);
     *dbp2 += fabs(ufinal[center_x])/2;

      /*find the probability of crossing the lower bound at this moment*/
     if ( lb_margin )
      { margin = lb_change ? lb_margin+lb_change[t] : lb_margin;
        for ( x=0,dbp2=Ptlow+r,*dbp2=0,dbp3=ufinal; (long)x < margin ; x++,dbp3++ )
         { *dbp2 += fabs(*dbp3);		*dbp3 = 0;	}		}

     	/*find the probability of crossing the upper bound at this moment*/
     if ( ub_margin )
      { margin = ub_change ? ub_margin-ub_change[t] : ub_margin;
        for ( x=0,dbp2=Pthigh+r,*dbp2=0,dbp3=ufinal+ulen-margin ; (long)x < margin ; x++,dbp3++ )
         { *dbp2 += fabs(*dbp3);		*dbp3 = 0;	}		}

     	/*now calculate the probability of not crossing the boundary upto this point in time*/
     for ( x=0,dbp2=Pt+r,*dbp2=0,dbp3=ufinal ; x < ulen ; x++,dbp3++ )
      *dbp2 += fabs(*dbp3);

      /*store the probability density at this moment in time in Pxt (the big matrix)*/
     if ( Pxt )
      { for ( x=0,dbp2=Pxt+ulen*t,dbp3=ufinal ; x < ulen ; x++,dbp2++,dbp3++ )
         *dbp2 = fabs(*dbp3);		}
   }

  		/*release the allocated memory*/
  if ( A )			mxFree ( A );
  if ( Wplus )		mxFree ( Wplus );
  if ( Wminus )	mxFree ( Wminus );
  if ( ubuf )		mxFree ( ubuf );

  return 1;
 }




/*******************
%
%
% Matlab call:
%
% 	[ufinal, Pt, Ptlow, Pthigh, Pg0, Pxt] = FP4 ( xmesh, uinit, k, sigma, [lb_change ub_change], [lb_margin; ub_margin], dt);
%
% Description:
%
% 	Propagates the probability density <uinit> on the spatial grid <xmesh> across time, using the diffusion
%    parameters k (drift rate) and sigma (diffusion standard deviation). The duration of the propagation is
%    defined by dt and size of [lb_change ub_change] (see below). Don't include Pxt in the output arguments
%    to keep memory management efficient and fast.
%
%
% input ->  [xmesh]						//a column vector specifying the spatial grid for probability density propagation in space
% 												//the values in xmesh are monotonically increasing and range from low values (lower bound - lb_margin + dx) to
%                             		//high values (upper bound + ub_margin - dx). xmesh should be regular. define it as:
%                                   //lower bound - lb_margin + dx : dx : upper bound + ub_margin - dx
%
%           [uinit]						//a column vector specifying the initial probabilities on the grid
%                                   //it is usually a delta function, such as:
%                                   //uinit = zeros(size(xmesh));   uinit(xmesh==0) = 1;
%
%           k								//a scalar specifying the drift rate
%
%           sigma							//a scalar specifying the standard deviation of the diffusion process
%
%           [lb_change ub_change]	//a matrix consisting of two columns, each column shows how the corresponding margin changes across time
%          									//the number of rows in this matrix defines the number of time steps for the propagation of the probability
%                                   //density. stepnum*dt is equal to the total propagation time
%                                   //make sure that the change in the bound height is valid: bound height should change so that the upper
%                                   //bound stay above zero and the lower bound stay below zero all the time, also neither of them can
%                                   //exceed the xmesh
%
%           [lb_margin; ub_margin]	//a column vector specifying the margin of the lower and upper boundaries
%          									//-in Fokker-Planck we normally get the probability of boundary crossing
%            								// but don't distinguish the probability of crossing each of the boundaries. To find the probability
%                         				// of crossing the upper and lower boundaries separately we should reduce xinit and increase xend by
%                                   // many folds of sigma and calculate how much of the probability falls beyond the actual target boundaries.
%                         				// lb_margin and ub_margin define the distance of the actual boundaries from xinit and xend (the starting
%                                   // and ending point of xmesh. When lb_margin and ub_margin are zero xinit and xend correspond to the boundaries.
%
%           [dt]							//the size of each time step, specifying the resolution of the temporal grid
%
%
% output->  [ufinal]						//a column vector, final probability density at the end of the propagation
%
% 				[Pt]							//a column vector, it is the survivor function: the probability of NOT crossing the bounds up to each moment
%                                   //during the propagation
%
%           [Ptlow, Pthigh]			//a matrix with two columns, the first column corresponds to the probability of crossing the lower bound
%          									//at each moment. the second column corresponds to the probability of crossing the upper bound at each moment
%
%           [Pg0]							//a column vector, the probability of staying above zero at each moment. it is useful when
%          									//the diffusion process terminates before bound crossing
%                                   //I have written the program so that it includes also half of the mass in the spatial grid location that
%                                   //corresponds to zero
%
%           [Pxt]							//a matrix with columns corresponding to the probability density at each moment during
%          									//the propagation. its calculation and storage requires a huge amount of memory. don't ask for
%                                   //it if you are not sure that your computer can handle it or not.
%
%
*******************/



void mexFunction ( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
 {
  double *xmesh, dx;
  double *uinit, *ufinal;
  unsigned long ulen;
  double *lb_change, *ub_change;
  double lb_margin, ub_margin;
  double k, sigma;
  double dt;
  unsigned long center_x;
  unsigned long stepnum;
  double *Pt, *Ptlow, *Pthigh, *Pg0, *Pxt=0;
  double *dbp1, *dbp2, dum;
  int i;
  long err;


			/*make sure all the necessary parameters are passed to the mex*/
  if ( nrhs != 7 )
   mexErrMsgTxt ( "Wrong number of arguments\n"\
                  " [ufinal, Pt, Pt_bound, Pg0, Pxt] = function (xmesh, uinit, k, sigma, bound_change, bound_margin, dt)" );


         /*get the spatial grid*/
  xmesh = mxGetPr(prhs[0]);
  			/*get the initial probability density*/
  uinit = mxGetPr(prhs[1]);
  			/*get the length of the spatial grid*/
  ulen = mxGetM(prhs[1]);
  			/*make sure that xmesh and uinit are both column vectors*/
  if ( mxGetN(prhs[0]) != 1 || mxGetN(prhs[1]) != 1 )
   mexErrMsgTxt ( "the spatial grid and the initial probability distribution should be column vectors" );
   		/*make sure that xmesh and uinit have the same size*/
  if ( ulen != mxGetM(prhs[0]) )
   mexErrMsgTxt ( "mismatch in the size of the spatial grid and the initial probability distribution" );
   		/*make sure that xmesh is arranged from the lower boundary to the upper boundary and not vice versa*/
  if ( xmesh[0] > xmesh[ulen-1] || xmesh[0] >= 0 || xmesh[ulen-1] <= 0 )
   mexErrMsgTxt ( "the spatial grid should be monotonically increasing\n\t\te.g., it can't go from the upper bound to the lower bound" );


   		/*get the dx and make sure that xmesh is a reqular grid, ie. dx is the same everywhere on the grid*/
  dx = xmesh[1]-xmesh[0];
  for ( i=1 ; i < ulen-1 ; i++ )
   if ( fabs((xmesh[i+1]-xmesh[i])-dx) > dx/100 )
    mexErrMsgTxt ( "the spatial grid should be regular\n\t\tdx should be constant everywhere on the grid" );


    		/*on the spatial grid find the center (where corresponds to zero)*/
  center_x = floor(-xmesh[0]/dx+0.5);


			/*get the drift rate*/
  k = mxGetScalar(prhs[2]);


  			/*get the diffusion standard deviation*/
  sigma = mxGetScalar(prhs[3]);


  			/*get the change in the position of the lower and upper boundaries*/
  dbp1 = mxGetPr(prhs[4]);
  if ( mxGetN(prhs[4]) != 2 )
   mexErrMsgTxt ( "bound_change should be a matrix with two columns, each of\n\t\twhich corresponding to the change of one of the boundaries across time" );
  stepnum = mxGetM(prhs[4]);
  lb_change = dbp1;
  ub_change = dbp1+stepnum;


  			/*get the margin for the lower and upper bounds*/
  dbp1 = mxGetPr(prhs[5]);
  lb_margin = dbp1[0];
  ub_margin = dbp1[1];
  			/*make sure that the size of margin vector is correct*/
  if ( mxGetNumberOfElements(prhs[5]) != 2 )
   mexErrMsgTxt ( "bound_margin should be a vector with two elements" );
   		/*make sure that bound changes do not eliminate the margins and also do not cause the upper bound go
           below zero or the lower bound go above zero*/
  for ( i=0,dbp1=lb_change,dbp2=ub_change ; i < stepnum ; i++,dbp1++,dbp2++ )
   if ( lb_margin+(*dbp1) < 3*sigma || ub_margin-(*dbp2) < 3*sigma ||
        lb_margin+xmesh[0]+(*dbp1)+dx >= xmesh[ulen-1]-ub_margin+(*dbp2) )
     mexErrMsgTxt ( "excessive bound change\n\t\tboundaries should not exceed the specified spatial grid,\n\t\tmargins should be larger than 3*sigma at all times\n\t\talso upper and lower bounds should not cross each other");


     		/*get the temporal resolution*/
  dt = mxGetScalar(prhs[6]);


  		/*summary of input parameters*/
/*
  mexPrintf("xinit:%f, xend:%f, dx:%f, stepnum: %d, dt: %f\n",xmesh[0],xmesh[ulen-1],dx,stepnum,dt);
*/

		/*allocate memory for the output parameters, check for memory allocation errors*/
  if ( (ufinal = (double *)mxCalloc(ulen,sizeof(double))) == NULL ||
       (Pt = (double *)mxCalloc(stepnum+1,sizeof(double))) == NULL ||
       (Ptlow = (double *)mxCalloc(stepnum+1,sizeof(double))) == NULL ||
       (Pthigh = (double *)mxCalloc(stepnum+1,sizeof(double))) == NULL ||
       (Pg0 = (double *)mxCalloc(stepnum+1,sizeof(double))) == NULL )
    mexErrMsgTxt ( "error in memory allocation" );
  if ( nlhs > 4 )
   if ( (Pxt = (double *)mxCalloc(ulen*stepnum,sizeof(double))) == NULL )
    mexErrMsgTxt ( "error in memory allocation" );


		/*propagate the initial probability density and get the final probability density, the survivor function
        and the probability of bound crossing at each moment*/
  err = !pde4 ( uinit, ulen, center_x, dx,
                lb_margin, ub_margin, lb_change, ub_change,
                stepnum, dt,
                k, sigma,
                ufinal, Pt, Ptlow, Pthigh, Pg0, Pxt );


		/*set the output arguments appropriately, this showed to be more tricky that I expected, the following code
        handles different number of outputs appropriately and the caller does not have to ask for all the possible
        outputs*/
  if ( nlhs > 0 )
   { for ( i = 0 ; i < ulen ; i++ )		ufinal[i] = fabs(ufinal[i]);
     plhs[0] = mxCreateDoubleMatrix(ulen,1,mxREAL);
     mxFree(mxGetPr(plhs[0]));
     mxSetPr(plhs[0],ufinal);		}
  else
   { mxFree(ufinal);					}

  if ( nlhs > 1 )
   { plhs[1] = mxCreateDoubleMatrix(stepnum+1,1,mxREAL);
     mxFree(mxGetPr(plhs[1]));
     mxSetPr(plhs[1],Pt);			}
  else
   { mxFree(Pt);						}

  if ( nlhs > 2 )
   { plhs[2] = mxCreateDoubleMatrix(stepnum+1,2,mxREAL);
     memcpy(mxGetPr(plhs[2]),Ptlow,(stepnum+1)*sizeof(double));
     memcpy(mxGetPr(plhs[2])+stepnum+1,Pthigh,(stepnum+1)*sizeof(double));	}
  mxFree(Ptlow);
  mxFree(Pthigh);

  if ( nlhs > 3 )
   { plhs[3] = mxCreateDoubleMatrix(stepnum+1,1,mxREAL);
     mxFree(mxGetPr(plhs[3]));
     mxSetPr(plhs[3],Pg0);			}
  else
   { mxFree(Pg0);						}

  if ( nlhs > 4 )
   { plhs[4] = mxCreateDoubleMatrix(ulen,stepnum,mxREAL);
     mxFree(mxGetPr(plhs[4]));
     mxSetPr(plhs[4],Pxt);			}
  else
   { if ( Pxt )	mxFree(Pxt);	}


  		/*generate an error message if the propagation of the probability density has failed*/
  if ( err )
   mexErrMsgTxt ( "error in solving the Fokker-Plnck equation" );
 }



