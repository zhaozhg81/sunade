package sunade.stat.general;

import cern.colt.list.DoubleArrayList;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;
import cern.jet.stat.Descriptive;

public class kernel_density_estimation{
    public DoubleArrayList x, y;
    public int n;
    public double h; 

    public kernel_density_estimation(DoubleArrayList ini_x)
    {
        // If the tuning parameter is not chosen, use the one provided by         
        n = ini_x.size();
        x=ini_x.copy();
        
        
        h = 1.06 * Math.sqrt(  Descriptive.sampleVariance(x, Descriptive.mean(x)) 
        )/Math.pow(n, 0.2);
        
    }
    
    public kernel_density_estimation(DoubleArrayList ini_x, double ini_h)
    {
        h = ini_h; 
        n = ini_x.size();
        x=ini_x.copy();
    }
    
    public void run_kde()
    {
        DRand rand_eng = new DRand(20);
        Normal norm_rv=new Normal(0, 1, rand_eng);
        
        y=x.copy();
        int i,j;
        for(i=0;i<n; i++)
        {
            y.set(i, 0);
            for(j=0;j<n;j++)
            {
               y.set(i, norm_rv.pdf( (x.get(i)-x.get(j) )/h )/h + y.get(i) );
            }
            y.set(i, y.get(i)/n);
        }
    }
    
        public void run_kde(DoubleArrayList new_x)
    {        
        DRand rand_eng = new DRand(20);
        Normal norm_rv=new Normal(0, 1, rand_eng);
        
        y=new_x.copy();
        int i,j;
        for(i=0;i<y.size(); i++)
        {
            y.set(i, 0);
            for(j=0;j<n;j++)
               y.set(i, norm_rv.pdf( (new_x.get(i)-x.get(j) )/h )/h + y.get(i) );
            y.set(i, y.get(i)/n);
        }
    }

}