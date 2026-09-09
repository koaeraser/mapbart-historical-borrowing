
class heterbart : public bart
{
  public:
   heterbart():bart() { }
   heterbart(size_t m):bart(m) { }
   void draw(double *sigma, std::vector<double>& allLsum,
             rn& gen, int shards=1);

   //----- Per-tree split of draw() so a calibration hook can fire
   // between heterdrmu() and the re-fit, exposing the partial residual
   // r[k] and the just-updated tree t[j].  Heteroscedastic counterpart
   // of bart::draw_per_tree in the Gaussian fork.  svec is a length-n
   // array of per-obs residual SDs (svec[k] = iw[k]*sigma in AFT).
   void draw_per_tree(size_t j, double *svec, rn& gen, int shards=1) {
      fit(t[j], xi, p, n, x, ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      heterbd(t[j], xi, di, pi, svec, nv, pv, aug, gen, shards);
      heterdrmu(t[j], xi, di, pi, svec, gen);
   }
   void finish_per_tree(size_t j) {
      fit(t[j], xi, p, n, x, ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }
};

//--------------------------------------------------
void heterbart::draw(double *sigma, 
                     std::vector<double>& allLsum,
                     rn& gen, int shards)
{
   size_t i=0;
   for(size_t j=0;j<m;j++) {
     // cout << "      Tree " << j << "      "<< endl;
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      if(heterbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen,shards)) i++;
      heterdrmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
      
      hetergetdiff(t[j],allLsum[j]);
   }
   accept=i;
}
