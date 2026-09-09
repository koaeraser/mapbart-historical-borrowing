#include <ctime>

class bart {
public:
   bart();
   bart(size_t m);

   void setdata(size_t p, size_t n, double *x, double *y, size_t nc=100);
   void setdata(size_t p, size_t n, double *x, double *y, int* nc);
   void setprior(double alpha, double beta, double tau)
      {pi.alpha=alpha; pi.mybeta = beta; pi.tau=tau;}

   xinfo& getxinfo() {return xi;}
   void setxinfo(xinfo& _xi);
   
   std::vector<size_t>& getnv() {return nv;}
   std::vector<double>& getpv() {return pv;}
   void setpv(double *varprob) {
     for(size_t j=0;j<p;j++) pv[j]=varprob[j];
   }
   tree& gettree(size_t i ) { return t[i];}
   
   double gettheta() {return theta;}

   void predict(size_t p, size_t n, double *x, double *fp);
   double f(size_t i) {return allfit[i];}
   double getaccept() {return accept;}

   //----- Calibration-hook accessors (read-only views into BART internals)
   // so a per-tree calibration callback can compute leaf-residual summaries
   // without touching protected members directly.  Used by cabart_ess.cpp.
   double  getr(size_t k)     const {return r[k];}
   double* getx_row(size_t k)       {return &x[k*p];}
   size_t  getp_()            const {return p;}
   size_t  getn_()            const {return n;}
   xinfo&  get_xi()                 {return xi;}

protected:
   size_t m; 
   std::vector<tree> t; 
   pinfo pi; 
   size_t p,n;
   double *x,*y;  
   xinfo xi; 
   double *allfit; 
   double *r;
   double *ftemp;
   dinfo di;
   bool dart,dartOn,aug,const_theta;
   double a,b,rho,theta,omega,accept; 
   double *grp;
   std::vector<size_t> nv;
   std::vector<double> pv, lpv;
};

//--------------------------------------------------
//constructor
bart::bart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false),accept(0.),grp(0) {}
bart::bart(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false),accept(0.),grp(0) {}

//--------------------------------------------------
void bart::setxinfo(xinfo& _xi)
{
   size_t p=_xi.size();
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     size_t nc=_xi[i].size();
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
   }
}
//--------------------------------------------------
void bart::setdata(size_t p, size_t n, double *x, double *y, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  this->setdata(p, n, x, y, nc);
  delete [] nc;
}

void bart::setdata(size_t p, size_t n, double *x, double *y, int *nc)
{
   this->p=p; this->n=n; this->x=x; this->y=y;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

   if(allfit) delete[] allfit;
   allfit = new double[n];
   predict(p,n,x,allfit);

   if(r) delete[] r;
   r = new double[n];

   if(ftemp) delete[] ftemp;
   ftemp = new double[n];

   di.n=n; di.p=p; di.x = &x[0]; di.y=r;
   for(size_t j=0;j<p;j++){
     nv.push_back(0);
     pv.push_back(1/(double)p);
   }
}
//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *fp)
{
   double *fptemp = new double[n];

   for(size_t j=0;j<n;j++) fp[j]=0.0;
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,fptemp);
      for(size_t k=0;k<n;k++) fp[k] += fptemp[k];
   }

   delete[] fptemp;
}

