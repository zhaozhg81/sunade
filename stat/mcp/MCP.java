/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author zhaozhg
 */

package sunade.stat.mcp;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.BooleanArrayList;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;
import java.util.Arrays;
import sunade.stat.general.QUICKSORT;
import sunade.stat.general.kernel_density_estimation;

public class MCP {
    
    // Public assessible variables.
    public BooleanArrayList decision;
    public String Method;
    public double alpha;
    
    // Private variables
    
    DoubleArrayList test_stat, p_value;

    int number_hypotheses, data_length;        
    int number_of_rejection;
    double threshold = 0;
 
    
    public int return_number_of_rejection()
    {
        return(number_of_rejection);
    }
    public void set_number_hypotheses(int nh){
        number_hypotheses=nh;
    }
    
    public MCP(DoubleArrayList data, String Is_P_value)
    {
        int i;
        if( "p_value".equals(Is_P_value))
        {
            number_hypotheses = data.size();
            p_value=data.copy();
            data_length=data.size();
            
        }else{
            number_hypotheses=data.size();
            test_stat=data.copy();
            data_length=data.size();
        }
        decision=new BooleanArrayList(data_length);
        for(i=0;i<data_length;i++)
            decision.add(false);

        alpha = 0.05;
        number_of_rejection=0;
        threshold=0;
        Method="None";
    }
    
    public MCP(DoubleArrayList data, String Is_P_value, double ini_alpha)
    {
        int i;
        if( "p_value".equals(Is_P_value))
        {
            data_length = data.size();
            p_value=data.copy();
            number_hypotheses=data_length;
        }else{
            data_length=data.size();
            test_stat=data.copy();
            number_hypotheses=data_length;
        }
        decision=new BooleanArrayList(data_length);
        for(i=0;i<data_length;i++)
            decision.add(false);

        alpha = ini_alpha;
        number_of_rejection=0;
        threshold=0;
        Method="None";
    }



    private void reset_decision()
    {
        int i;
        for(i=0; i<data_length; i++)        
            decision.set(i, false);
        number_of_rejection=0;
    }
    
    private void run_bh()
    {
        int i;
        DoubleArrayList p_value_sort;
        p_value_sort=p_value.copy();
        p_value_sort.sort();
        for(i=data_length; i>0; i--)
        {
            if( p_value_sort.get(i-1) <= (alpha*i)/number_hypotheses )
            {
                number_of_rejection=i;
                threshold = p_value_sort.get(i-1);
                break;
            }
        }
        for(i=0; i<data_length; i++)
            decision.set(i, ( p_value.get(i) <= threshold ) );
    }
    
        private void run_naive()
    {
        int i;
        number_of_rejection=0;
        for(i=0; i< data_length; i++ )
        {
            if( p_value.get(i) < alpha )
            {
                number_of_rejection= number_of_rejection+1;
                decision.set(i,true);
            }else{
                decision.set(i,false);
            }
        }
    }
    
        private void run_bonferroni()
    {
        int i;
        number_of_rejection=0;
        for(i=0; i< data_length; i++ )
        {
            if( p_value.get(i) < alpha/number_hypotheses )
            {
                number_of_rejection= number_of_rejection+1;
                decision.set(i,true);
            }else{
                decision.set(i,false);
            }
        }
    }

    
        private void run_by()
    {
        int i;
        DoubleArrayList p_value_sort;
        p_value_sort=p_value.copy();
        double total=0;
        for(i=1; i<= data_length; i++)
            total = total + Math.pow(i, -1);

        p_value_sort.sort();
        for(i=data_length; i>0; i--)
        {
            if( p_value_sort.get(i-1) <= (alpha*i/number_hypotheses/total) )
            {
                number_of_rejection=i;
                threshold = p_value_sort.get(i-1);
                break;
            }
        }
        for(i=0; i<data_length; i++)
            decision.set(i, ( p_value.get(i) <= threshold ) );
    }

        private void run_locfdr()
        {
            DRand rand_eng = new DRand(20);
            Normal norm_rv=new Normal(0, 1, rand_eng);

            int i;
            DoubleArrayList esti_marginal, loc_fdr, sort_loc_fdr;
            loc_fdr = test_stat.copy();
            kernel_density_estimation kde;
            kde = new kernel_density_estimation(test_stat.copy(), 1);
            kde.run_kde();
            esti_marginal = kde.y.copy();
            for(i=0;i<data_length; i++)
                loc_fdr.set(i, norm_rv.pdf( test_stat.get(i) )/esti_marginal.get(i) );
        
            sort_loc_fdr = loc_fdr.copy();
            sort_loc_fdr.sort();
            double cum_sum=0;
            for(i=0; i<data_length; i++)
            {
                cum_sum = cum_sum + sort_loc_fdr.get(i);
                if( cum_sum > (i+1)*alpha ){
                    number_of_rejection=i;
                    if( i >0 ){
                        threshold = sort_loc_fdr.get(i);
                    }else{
                        threshold =0;
                    }                        
                    break;
                }
            }
            for(i=0; i<data_length; i++)
                decision.set(i, loc_fdr.get(i) <= threshold );
        }
    
    public void run_analysis(String user_method) throws sunade_MCP_Exception
    {
        switch( user_method.toLowerCase() )
        {
            case "bh":
                reset_decision();
                run_bh();
                break;
        
            case "by":                
                reset_decision();
                run_by();
                break;
                
            case "locfdr":
                reset_decision();
                run_locfdr();
                break;
                
            case "naive":
                reset_decision();
                run_naive();
                break;
                
            case "bonferroni":
                reset_decision();
                run_bonferroni();
                break;
                
            default:
                throw new sunade_MCP_Exception("Wrongly specified name for the procedure.");
        }        
    }    
    
    public static boolean[] cut_BH(double[] pvalues, double alpha)
    {
        int p = pvalues.length;
        boolean[] decision = new boolean[p];
        return(decision);
    }
    
    public static boolean[] cut_locfdr(double[] lfdr, double alpha)
    {
        int i;
        int p= lfdr.length;
        boolean[] decision = new boolean[p];
               
        double [] sorted=new double[p];
        sorted = Arrays.copyOf(lfdr, p);
        QUICKSORT.quicksort(sorted);
        double cum_sum =0;
        int number_rej=0;
        
        for(i=0; i<p; i++)
        {
            cum_sum = cum_sum + sorted[i];
            if( cum_sum/(i+1) > alpha)
                break;
            number_rej = number_rej + 1;
        }
        double thresh=0;
        if(number_rej==0)
        {
            for(i=0;i<p;i++)
                decision[i]=false;
        }else{
            thresh = sorted[number_rej-1];
            for(i=0;i<p;i++)
                decision[i] = lfdr[i] <= thresh;
        }
        
        return(decision);
    }
    
    
    public static boolean[] cut_locfdr(double[] lfdr, double alpha, int[] Rg)
    {
        // Here the Rg argument plays as a weight parameter.
        // This functions is used in GATE 2.
        int i;
        int p= lfdr.length;
        boolean[] decision = new boolean[p];
               
        double [] sorted=new double[p];
        int [] Rg_sorted = new int[p];
        sorted = Arrays.copyOf(lfdr, p);
        Rg_sorted = Arrays.copyOf(Rg, p);
        
        QUICKSORT.quicksort(sorted, Rg_sorted);
        double cum_sum =0, cum_sum_Rg=0;
        int number_rej=0;
        
        for(i=0; i<p; i++)
        {
            cum_sum = cum_sum + sorted[i]*Rg_sorted[i];
            cum_sum_Rg = cum_sum_Rg + Rg_sorted[i];
            if( cum_sum/cum_sum_Rg > alpha)
                break;
            number_rej = number_rej + 1;
        }
        double thresh=0;
        if(number_rej==0)
        {
            for(i=0;i<p;i++)
                decision[i]=false;
        }else{
            thresh = sorted[number_rej-1];
            for(i=0;i<p;i++)
                decision[i] = lfdr[i] <= thresh;
        }
        
        return(decision);
    }

            
}
