/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat.mcp;

import java.util.ArrayList;

/**
 *
 * @author zhaozhg
 */
public class GATE_Data {
    int G;
    int[] group_sizes;
    int total_hypotheses;
    ArrayList<double []> data = new ArrayList<>();
    
    public GATE_Data(int ini_G, int[] ini_group_sizes, double[] ini_all_data)
    {
        int i,j;
        int cur_ind=0;
        G = ini_G;
        group_sizes = new int [G];
        
        double[] data_group_i;

        for(i=0;i<G; i++)
        {
            group_sizes[i] = ini_group_sizes[i];
            data_group_i = new double[group_sizes[i]];
            for(j=0; j<group_sizes[i]; j++)
                data_group_i[j] = ini_all_data[cur_ind + j];
            data.add(data_group_i);
            cur_ind = cur_ind + group_sizes[i];
        }
        total_hypotheses = ini_all_data.length;
    }
    public int get_G(){ return(G);}
    public int[] get_group_sizes(){ return(group_sizes);}
    public double[] get_data_given_group(int g){ return( data.get(g) ); }
    public int get_total_hypotheses(){ return(total_hypotheses); }
}
