
bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, 
	     std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen, int shards)
{
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);

      // cout << "v = " << v << endl;
      // cout << "c = " << c << endl;
      //--------------------------------------------------
      //compute sufficient statistics
      size_t nr,nl; //counts in proposed bots
      double bl,br; //sums of weights
      double Ml, Mr; //weighted sum of y in proposed bots
      hetergetsuff(x,nx,v,c,xi,di,
                   nl,bl,Ml,
                   nr,br,Mr,
                   sigma,shards);
      // cout << "nl = " << nl << endl;
      // cout << "nr = " << nr << endl;
      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
         lhl = heterlh(bl,Ml,pi.tau);
         lhr = heterlh(br,Mr,pi.tau);
         lht = heterlh(bl+br,Ml+Mr,pi.tau);
   
         alpha=1.0;
         lalpha = log(pr) + (lhl+lhr-lht); 
         lalpha = std::min(0.0,lalpha);
      }
      // cout << "lhl = " << lhl << endl;
      // cout << "lhr = " << lhr << endl;
      // cout << "lht = " << lht << endl;
      //--------------------------------------------------
      //try metrop
      double mul,mur; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
        // mul = heterdrawnodemu(bl,Ml,pi.tau,gen);
        // mur = heterdrawnodemu(br,Mr,pi.tau,gen);
        mul = gen.normal();
        double mul_ = gen.normal();
        mur = gen.normal();
        double mur_ = gen.normal();
        
         x.birthp(nx,v,c,mul,mur);
	 nv[v]++;
         return true;
      } else {
         return false;
      }
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree::tree_p nx; //nog node to death at
      dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      double br,bl; //sums of weights
      double Ml, Mr; //weighted sums of y
      hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigma, shards);

      //--------------------------------------------------
      //compute alpha
      double lhl, lhr, lht;
      lhl = heterlh(bl,Ml,pi.tau);
      lhr = heterlh(br,Mr,pi.tau);
      lht = heterlh(bl+br,Ml+Mr,pi.tau);

      double lalpha = log(pr) + (lht - lhl - lhr);
      lalpha = std::min(0.0,lalpha);

      //--------------------------------------------------
      //try metrop
      //double a,b,s2,yb;
      double mu;
      if(log(gen.uniform()) < lalpha) {
         // mu = heterdrawnodemu(bl+br,Ml+Mr,pi.tau,gen);
         mu = gen.normal();
        double mu_ = gen.normal();
	 nv[nx->getv()]--;
         x.deathp(nx,mu);
         return true;
      } else {
         return false;
      }
   }
}

