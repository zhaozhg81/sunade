/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat.general;

import sunade.stat.bet.BET_DataFile;
import sunade.stat.bet.BET_MultiVariate;
import sunade.stat.bet.BET_Pairwise;

/**
 *
 * @author zhaozhg
 */
public class R_Interface {
    public static double [] multi_bet_R(double [] R_Us, int dim, int sample_size, int depth)
    {
        int i,j;
        double[] result= new double[4+ dim];        
        double[][] Us = new double[sample_size][dim];
        for(i=0;i<sample_size;i++)
            for(j=0;j<dim;j++)
                Us[i][j] = R_Us[ i + j*sample_size];
        
        
        // This part is to debug if the data inputted from R is correct.
        /*
        System.out.println("dim="+dim+"\n");
        System.out.println("sample_size="+sample_size+"\n");
        System.out.println("depth="+depth+"\n");
        for(i=0;i<3;i++)
        {
            for(j=0;j<dim;j++)
            {
                System.out.printf("data[%d,%d]=%f,", i,j,Us[i][j]);
            }
            System.out.println("\n");
        }
        */
        
        BET_MultiVariate bet_high = new BET_MultiVariate(dim, depth, sample_size, Us);
        bet_high.get_As_Bs();
        bet_high.get_result();
        
        result[0] = bet_high.return_Z_stat();
        result[1] = bet_high.return_bonf_p_value();
        result[2] = bet_high.return_p_value();
        result[3] = (double) bet_high.return_Max_Diff();
        for(i=0;i<dim;i++)
            result[4+i] = bet_high.return_MAX_A_ind(i);
        
        return(result);
    }

}
