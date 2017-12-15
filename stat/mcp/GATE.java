/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat.mcp;

import cern.colt.list.DoubleArrayList;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class GATE {    
    int DEBUG=1;
//    GATE_Data_each_group[] gate_data = new GATE_Data_each_group[10];
    ArrayList<GATE_each_group> all_groups= new ArrayList<>(); //
    int G; // The total number of groups.
    double alpha; // The alpha level for the GATE
    int L=2; //the component for the Gaussian mixture
    GATE_Data data;    // the object to store all the data
    double pi0, pi1, pi20, pi21; // The parameters corresponding to all the proportions.
    double[] cL, muL, sigmaL;
    int[] group_sizes;
    int total_hypotheses;
    int rejected_group;
    int total_rejection;
    ArrayList<int[]> rejection_indicies = new ArrayList<>();
    int total_rej, total_group_rej;
    int[] number_rej_each_group;

    // Private variables
    
    DoubleArrayList[] test_stat, p_value;

    int number_hypotheses, data_length;        
    int number_of_rejection;
    double threshold = 0;

    public GATE(GATE_Data ini_data){
        data=ini_data;
        G = data.get_G();
        group_sizes = new int[G];
        group_sizes = data.get_group_sizes();
        int i;
        total_hypotheses=0;
        for(i=0; i<G; i++)
        {
            GATE_each_group gate= new GATE_each_group( data.get_data_given_group(i), L );
            all_groups.add(gate);        
            total_hypotheses = total_hypotheses + group_sizes[i];
        }
        cL = new double[L]; muL = new double[L]; sigmaL = new double[L];
        pi0=0.5; pi1=1-pi0;
        pi20=0.4; pi21=1-pi20;
        alpha = 0.05;
        cL[0]=0.3; cL[1]=0.7;
        muL[0]=1; muL[1]=-1;
        sigmaL[0]=Math.sqrt(1); sigmaL[1]=Math.sqrt(1);
    }

    public void calculate_Lfdr_all_groups()
    {
        int i;
        for(i=0; i<G; i++)
        {
            all_groups.get(i).set_pis(pi1, pi21);
            all_groups.get(i).set_mixed_gaussian_parameter(muL, sigmaL, cL);
            all_groups.get(i).cal_locfdr();
        }
    }
    
    public void calculate_Lfdr_all_groups(double pi1_tmp, double pi21_tmp, double[] cL_tmp, double[] muL_tmp, double[] sigmaL_tmp)
    {
        int i;
        for(i=0; i<G; i++)
        {
            all_groups.get(i).set_pis(pi1_tmp, pi21_tmp);
            all_groups.get(i).set_mixed_gaussian_parameter(muL_tmp, sigmaL_tmp, cL_tmp);
            all_groups.get(i).cal_locfdr();
        }
    }

    
    public void em(double DELTA, boolean sigma_KNOWN)
    {
        double diff = 1;
        double pi0_old, pi1_old, pi20_old, pi21_old; // The parameters corresponding to all the proportions.
        double[] cL_old, muL_old, sigmaL_old;
        double pi0_new, pi1_new, pi20_new, pi21_new; // The parameters corresponding to all the proportions.
        double[] cL_new, muL_new, sigmaL_new;
        double[] all_lfdr_g = new double[G];
        int i,j,l;
        
        pi0_old = pi0; pi1_old =pi1; pi20_old =pi20; pi21_old=pi21; 
        cL_old = new double[L]; muL_old= new double[L]; sigmaL_old = new double[L];
        cL_new = new double[L]; muL_new= new double[L]; sigmaL_new = new double[L];
        for(i=0; i<L; i++) { cL_old[i] = cL[i]; muL_old[i]=muL[i]; sigmaL_old[i]=sigmaL[i];}
        
        while( diff> DELTA )
        {
            pi1_new = 0; pi0_new = 1-pi1_new; pi21_new =0; pi20_new=1-pi21_new;
            calculate_Lfdr_all_groups(pi1_old, pi21_old, cL_old, muL_old, sigmaL_old);
            for(i=0; i<G; i++) all_lfdr_g[i]=all_groups.get(i).get_Lfdr_g();

            for(i=0;i<L;i++) { cL_new[i]=0; muL_new[i]=0; sigmaL_new[i]=0;}
            // The following code update pi1_new
            for(i=0;i<G; i++)
                pi1_new = pi1_new + ( 1- all_lfdr_g[i] );
            pi1_new = pi1_new /G; pi0_new = 1-pi1_new;
            
            // The following code update cl_new
            double den=0, num=0, tmp=0;
            for(l=0;l<L; l++)
            {
                den=0; num=0; tmp=0;
                for(i=0; i<G; i++)
                {            
                    for(j=0; j<all_groups.get(i).get_group_size(); j++)
                    {
                        num = num + all_groups.get(i).get_prob_m_il()[j][l];
                        den = den + ( 1- all_lfdr_g[i] )* (1-all_groups.get(i).get_Lfdr_gj(j) );
                        tmp = tmp + all_groups.get(i).get_prob_m_il()[j][0] + all_groups.get(i).get_prob_m_il()[j][1];

                    }
                }
                cL_new[l] = num/den;
            }
            
            // The following code update pi21_new;
            double A=0, B=0;
            for(i=0; i<G; i++)
            {
                for(j=0; j<all_groups.get(i).get_group_size(); j++)
                {
                    A = A + ( 1-all_lfdr_g[i]) * all_groups.get(i).get_Lfdr_gj(j);
                    B = B + ( 1-all_lfdr_g[i]) * (1-all_groups.get(i).get_Lfdr_gj(j));
                }
            }
            pi21_new = solve_pi21_new(A, B, all_lfdr_g);
            // The following code update both mul and sigmaL;
            for(l=0;l<L; l++)
            {
                num=0; den=0;
                for(i=0;i<G; i++)
                {
                    for(j=0;j<all_groups.get(i).get_group_size(); j++)
                    {
                        num = num + all_groups.get(i).get_data(j) * all_groups.get(i).get_prob_m_il()[j][l];
                        den = den + all_groups.get(i).get_prob_m_il()[j][l];
                    }
                }
                muL_new[l] = num/den;
            }
            // update sigmaL
            if( sigma_KNOWN==false){
                for (l = 0; l < L; l++) {
                    num = 0;
                    den = 0;
                    for (i = 0; i < G; i++) {
                        for (j = 0; j < all_groups.get(i).get_group_size(); j++) {
                            num = num + Math.pow(all_groups.get(i).get_data(j) - muL_new[l], 2) * all_groups.get(i).get_prob_m_il()[j][l];
                            den = den + all_groups.get(i).get_prob_m_il()[j][l];
                        }
                    }
                    sigmaL_new[l] = Math.sqrt(num / den);
                }
            }else{
                for(l=0; l<L; l++)
                    sigmaL_new[l] = sigmaL_old[l];
            }
            diff=0;
            diff= Math.pow(pi21_new-pi21_old, 2) + Math.pow( pi1_new - pi1_old, 2);
            for(l=0; l<L ; l++)
            {
                diff= diff + Math.pow(cL_new[l]-cL_old[l], 2);
                diff= diff + Math.pow(muL_new[l]-muL_old[l], 2);
                diff= diff + Math.pow(sigmaL_new[l]-sigmaL_old[l], 2);   
                pi1_old = pi1_new; pi0_old = pi0_new; pi20_old =pi20_new; pi21_old =pi21_new;
                
                cL_old[l] = cL_new[l];
                muL_old[l] = muL_new[l];
                sigmaL_old[l] = sigmaL_new[l];
            }
            if(DEBUG==1)
                System.out.println("delta=" + Double.toString(diff) +"\n");
        }
        // Finish the loop and set the value for all the parameters.
        pi1 = pi1_old; pi0 = 1-pi1; pi21 = pi21_old; pi20=1-pi21;
        for(l=0;l<L; l++)
        { cL[l] = cL_old[l]; muL[l] = muL_old[l]; sigmaL[l] = sigmaL_old[l];}
        
    }
    
    private double solve_pi21_new(double A, double B, double[] all_lfdr_g)
    {
        int increment=1000;
        double[] x0 = new double[increment-100];
        double[] y0 = new double[increment-100];
        double x;
        int i,j;
        double MAX;
        int max_ind;
        
        for(i=0;i<(increment-100);i++)
        {
            x0[i] = (double) (i+1)/(increment+1) + 0.1;
            x=x0[i];
            
            y0[i] = A * Math.log(1-x) + B *Math.log(x);
            for(j=0; j<G; j++)
            {
                y0[i] = y0[i] - Math.log( 1- Math.pow(1-x, group_sizes[j])) * group_sizes[j] * (1 - all_lfdr_g[j] );
            }
        }
        MAX=y0[0]; max_ind=0;
        for(i=1;i<(increment-100);i++)
        {
            if(y0[i] > MAX)
            {
                MAX=y0[i]; max_ind = i;
            }
        }
        return(x0[max_ind]);
    }
    
    public void perform_gate_1()
    {
        int i,j;
        int cur_ind=0;
        calculate_Lfdr_all_groups();
        double[] all_locfdr = new double[ total_hypotheses ];
        boolean[] all_decision = new boolean[ total_hypotheses ];
        
        // Put all the local fdr scores together
        for(i=0; i<G; i++)
        {
            for(j=0; j<group_sizes[i]; j++)
            {
                all_locfdr[ cur_ind ]= 1- (1-all_groups.get(i).get_Lfdr_g() )*( 1- all_groups.get(i).get_Lfdr_gj(j));
                cur_ind = cur_ind +1 ;
            }
        }
        // Apply the local fdr approach on these local fdr scores.
        
        all_decision = MCP.cut_locfdr(all_locfdr, alpha);
        rejected_group=0;
        total_rejection = 0;
        
        cur_ind=0;
        boolean group_dec;
        int[] temp= new int[2];
        // Set all the rules for all the groups.
        for(i=0;i<G;i++)
        {
            group_dec=false;
            for(j=0;j<group_sizes[i]; j++)
            {
                all_groups.get(i).set_wg_decision(j, all_decision[cur_ind]);
                if( all_decision[cur_ind]== true)
                    { 
                        group_dec=true; total_rejection +=1;
                        temp[0] = i+1;
                        temp[1] = j+1;
                        rejection_indicies.add(temp);
                    }
                cur_ind = cur_ind + 1;
            }
            if( group_dec== true)
                rejected_group += 1;
            all_groups.get(i).set_bg_decision(group_dec);
        }
        summarize_result();
        System.out.println("Done the calculation.");
    }
    
    public ArrayList<int[]> get_rejection_indicies(){ return(rejection_indicies);}
    
    public void perform_gate_2(double eta)
    {
        int i,j, mg;
        int cur_ind=0;
        calculate_Lfdr_all_groups();
        boolean[] wg_decision;
        double[] cum_group_lfdr = new double[ G ];
        int[] Rg = new int[G];
        
        // Step 1. Calculate the "potential" within group decision.
        // Calculte the cumulage average of local fdr scores within the groups.
        for(i=0; i<G; i++)
        {
            mg = all_groups.get(i).get_group_size();
            wg_decision = new boolean[mg];
            wg_decision = MCP.cut_locfdr( all_groups.get(i).get_Lfdr_gj_all(), eta );
            all_groups.get(i).set_wg_decision_all(wg_decision);
            Rg[i] = 0;
            cum_group_lfdr[i] =0;
            for(j=0; j<mg; j++)
            {
                if(  wg_decision[j]== true )
                {
                    Rg[i] = Rg[i] + 1;
                    cum_group_lfdr[i] = cum_group_lfdr[i] + all_groups.get(i).get_Lfdr_gj(j);
                }
            }
            if ( Rg[i] >0 ){
                cum_group_lfdr[i] = cum_group_lfdr[i]/Rg[i];
            }
        }
        // Step 2. Determine the between-group decision.
        double[] fdr_g_star = new double[G];
        boolean[] group_dec = new boolean[G];
        for(i=0; i<G; i++ )
        {
            fdr_g_star[i] = 1 - (1-cum_group_lfdr[i]) * (1-all_groups.get(i).get_Lfdr_g());
        }
        group_dec = MCP.cut_locfdr(fdr_g_star, alpha, Rg);
        
        // Step 3. Combine the decision on the group-level and potentially rejected hypotheses
        //         in each group. We come to the final rejections of all the hypotheses.
        
        for(i=0; i<G; i++)
        {
            if( group_dec[i] == true ){
                all_groups.get(i).set_bg_decision(true);
            }else{
                all_groups.get(i).set_bg_decision(false);
                for(j=0; j< group_sizes[i]; j++)
                {
                    all_groups.get(i).set_wg_decision(j, false);
                }
            }
        }
        
        summarize_result();
        System.out.println("Done the calculation of GATE 2.");        
    }
            
            
    public void perform_gate_3()
    {
        
    }
    
    public void perform_gate_4(double eta, int increment)
    {
        int i,j;
        double[] alpha_prime = new double[increment];
        double[] PFDR_Sel = new double[increment];
        double [] Lfdr_gj;
        double Lfdr_g;
        boolean[] decision;
        
        
        calculate_Lfdr_all_groups();
        for(i=0;i<increment;i++)
            alpha_prime[i] = alpha*(i+1)/increment;
        // Step one: select the groups such that the cumulative average of the local fdr
        //      scores on the group level is controlled at eta level.        
        double [] all_lfdr_g = new double[G];
        boolean [] selection_index = new boolean[G];
        for(i=0; i<G; i++)
            all_lfdr_g[i] = all_groups.get(i).get_Lfdr_g();
        selection_index = MCP.cut_locfdr(all_lfdr_g, eta);
        
        // Step 2. Find the maximum alpha' such that the selective PFDR is controlled at alpha level.
        double alpha_star=0;
        for(i=0;i<increment;i++)
        {
            alpha_prime[i] = alpha*(i+1)/increment;
            PFDR_Sel[i] = calculate_PFDR_Sel( alpha_prime[i], selection_index);
            if( PFDR_Sel[i] <= alpha )
                alpha_star = alpha_prime[i];
        }
        if( (alpha_star ==0)| (count_number_true(selection_index) == 0) )
        {
            for(i=0; i<G; i++)
            {
                all_groups.get(i).set_bg_decision(false);
                for(j=0; j<all_groups.get(i).get_group_size(); j++)
                    all_groups.get(i).set_wg_decision(j, false);
            }
        }else{
            for(i=0;i<G; i++)
            {
                if (selection_index[i] == true) {
                    Lfdr_gj = new double[group_sizes[i]];
                    decision = new boolean[group_sizes[i]];
                    Lfdr_gj = all_groups.get(i).get_Lfdr_gj_all();                    
                    decision = MCP.cut_locfdr(Lfdr_gj, alpha);
                    if( count_number_true(decision)==0 )
                    {
                        all_groups.get(i).set_bg_decision(false);
                        for(j=0; j<all_groups.get(i).get_group_size(); j++)
                            all_groups.get(i).set_wg_decision(j, false);
                    }else{
                        all_groups.get(i).set_bg_decision(true);
                        for(j=0; j<all_groups.get(i).get_group_size(); j++)
                            all_groups.get(i).set_wg_decision(j, decision[j]);
                    }
                }else{
                    all_groups.get(i).set_bg_decision(false);
                    for(j=0; j<all_groups.get(i).get_group_size(); j++)
                        all_groups.get(i).set_wg_decision(j, false);
                }
            }        
        }
        summarize_result();
    }
    
    private void summarize_result()
    {
        int i,j;
        int number_rej;
        
        number_rej_each_group = new int[G];
        boolean[] temp;
                
        number_rej = 0; total_rej=0; total_group_rej=0;
        for (i = 0; i < G; i++) {
            if ( all_groups.get(i).get_bg_decision() == true ) {
                total_group_rej +=1;
                temp = new boolean[group_sizes[i]];
                temp = all_groups.get(i).get_wg_decision();
                for (j = 0; j < group_sizes[i]; j++) {
                    if (temp[j] == true) {
                        number_rej_each_group[i] += 1;
                    }
                }
                total_rej = total_rej + number_rej_each_group[i];
            } else {
                number_rej_each_group[i] = 0;
            }
        }

    }
    
    public void save_result(String filename) throws IOException
    {
        
        try {
            FileWriter output = new FileWriter(filename, true);
            
            int i,j;
            String to_be_write;
            for(i=0;i<G;i++)
            {
                for(j=0; j< all_groups.get(i).get_group_size(); j++)
                {
                    to_be_write = String.format("%d,%d,%.10f,%b,%b\n", i+1, j+1,
                            all_groups.get(i).get_data(j), all_groups.get(i).get_bg_decision(),
                            all_groups.get(i).get_wg_decision()[j]
                    );
                    output.write( to_be_write );
                }
            }
            output.close();
        }catch (IOException ex) {
            System.err.println(ex);
        }
       
    }
    
    
    private double calculate_PFDR_Sel(double alpha, boolean[] selection_index)
    {
        int i;
        double [] Lfdr_gj;
        double Lfdr_g;
        boolean[] decision;
        double PFDR_Sel = 0;
        int Rg;
        int Group_sel = 0;

        if (count_number_true(selection_index) > 0) {
            for (i = 0; i < G; i++) {
                if (selection_index[i] == true) {
                    Lfdr_gj = new double[group_sizes[i]];
                    decision = new boolean[group_sizes[i]];
                    Lfdr_gj = all_groups.get(i).get_Lfdr_gj_all();
                    Lfdr_g = all_groups.get(i).get_Lfdr_g();
                    decision = MCP.cut_locfdr(Lfdr_gj, alpha);
                    Rg = count_number_true(decision);
                    if (Rg > 0) {
                        Group_sel = Group_sel + 1;
                        // PFDR_Sel = PFDR_Sel + (1 - Lfdr_g) * (1 - sum_of_weighted_array(Lfdr_gj, decision) / Rg);
                        PFDR_Sel = PFDR_Sel + Lfdr_g + ( 1-Lfdr_g) * sum_of_weighted_array(Lfdr_gj, decision)/Rg;
                    }
                }
            }
            PFDR_Sel = PFDR_Sel/Group_sel;
        } else {
            PFDR_Sel = 1;
        }
        return (PFDR_Sel);
    }
    
    private int count_number_true( boolean[] dec)
    {
        int i;
        int R=0;
        for(i=0; i<dec.length; i++)
            if( dec[i] == true )
                R = R + 1;
        return(R);
    }
    private double sum_of_array( double[] array)
    {
        double sum=0;
        int i;
        for(i=0; i<array.length; i++)
            sum=sum+ array[i];
        return(sum);
    }
    private double sum_of_weighted_array( double[] array, double[] weight)
    {
        double sum=0;
        int i;
        for(i=0; i<array.length; i++)
            sum=sum+ array[i]*weight[i];
        return(sum);
    }

    private double sum_of_weighted_array( double[] array, boolean[] weight)
    {
        double sum=0;
        int i;
        for (i = 0; i < array.length; i++) {
            if (weight[i] == true) {
                sum = sum + array[i];
            }
        }
        return(sum);
    }
    
    
    public boolean get_bg_rejection(int g){  return( all_groups.get(g).get_bg_decision() ); }
    public boolean [] get_wg_rejection(int g){ return( all_groups.get(g).get_wg_decision() ); }
    public int get_total_rej(){ return(total_rej);}
    public int get_group_rej(){ return(total_group_rej);}
    public int[] get_rej_groups(){ return( number_rej_each_group); }
    public int get_total_group() { return(G); }

    public void set_number_hypotheses(int nh){
        number_hypotheses=nh;
    }
    
 
    public static void main(String[] args) throws IOException
    {
        int G,i,j;
        int[] group_sizes;
        double[] data_all;
        
        String current = new java.io.File( "." ).getCanonicalPath();
        System.out.println("Current dir:"+current);
        
        String datafile = "R/AYP.csv";
        GATE_ReadFile gate_read = new GATE_ReadFile(datafile);
        G = gate_read.get_G();
        group_sizes = gate_read.get_group_sizes();
        data_all = gate_read.get_data_all();
        
        GATE_Data gate_data = new GATE_Data(G, group_sizes, data_all);
        GATE gate = new GATE(gate_data);
        gate.em(0.001, true);
        // gate.perform_gate_1();
        //gate.perform_gate_2(0.05);
        gate.perform_gate_4(0.01,100);

        int total_rej;
        int[] number_rej_each_group = new int[G];

        total_rej=gate.get_total_rej();
        number_rej_each_group = gate.get_rej_groups();


        System.err.println("Done the calculation.");

    }
    
}
