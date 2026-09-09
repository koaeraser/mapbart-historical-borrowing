
class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;q=0;}
   size_t p; 
   size_t n;  
   double *x; 
   double *y; 
   int *q;
};

class pinfo
{
public:
   pinfo(): pbd(1.0),pb(.5),alpha(.95),mybeta(2.0),tau(1.0) {}
   double pbd; 
   double pb;  
   double alpha;
   double mybeta;
   double tau;
};
