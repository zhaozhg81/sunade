/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat.mcp;

import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;

/**
 *
 * @author zhaozhg
 */

public class GATE_each_group{    
    int group_size;
    double[] test_stat;// the value of test statistic
    double[] Lfdr_gj; //P(\theta_{ij}=0｜data, \theta_i = 1)
    double Lfdr_g; // P(\theta_i=0|data)
    boolean[] wg_decision; //Within-group decision
    boolean bg_decision; // Between-group decision
    double pi1, pi21, pi0, pi20; // pi1, pi0 are the proportion on the group level.
    // pi1 is the proportion of significant groups.
    // pi21 is the proportion of a hypotheses to be significant within a significant group.
    int L; // Number of components in the Gaussian mixture
    double[] muL, sigmaL, cL;// Parameters corresponding to the Gaussian mixture
    double[] f0x, f1x, fx;
    double[][] prob_m_il; // The variable stores the components for the Gaussian mixture.
                          // It corresponding to the following quantities in Liu, Sarkar and Zhao(2016).
                          // P(theta_g=1, theta_{i|g}=1,m_{i|g}=l|data)
    
    
    public double get_Lfdr_g(){ return(Lfdr_g);}
    public double get_Lfdr_gj(int i){ return( Lfdr_gj[i]); }
    public boolean get_bg_decision(){ return(bg_decision);}
    public int get_group_size(){ return(group_size);}
    public double[] get_Lfdr_gj_all(){ return(Lfdr_gj); }
    public double get_data(int i){ return( test_stat[i]); }    
    public boolean[] get_wg_decision(){ return( wg_decision);}
    public double[][] get_prob_m_il() { return( prob_m_il); }
    
    
    public void set_bg_decision(boolean dec){ bg_decision=dec;}
    public void set_wg_decision(int i, boolean dec){ wg_decision[i]=dec;}
    public void set_wg_decision_all( boolean[] dec){ wg_decision=dec; }
    
    public GATE_each_group(double[] ts, int ini_L)
    {
        int i;
        group_size = ts.length;
        test_stat = new double [group_size];
        Lfdr_gj = new double[group_size];
        wg_decision = new boolean[group_size];
        f0x = new double[group_size];
        f1x = new double[group_size];
        fx = new double[group_size];
        for(i=0; i < group_size; i++)
        {
            test_stat[i] = ts[i];
            Lfdr_gj[i] = 0;
            wg_decision[i] = false;
        }
        bg_decision = false;        
        L = ini_L;
        prob_m_il = new double[group_size][L];
    }
    public void set_pis(double p1, double p21)
    {
        pi1 = p1;
        pi21= p21;
        pi0=1-pi1;
        pi20=1-pi21;
    }

    public void set_mixed_gaussian_parameter(double[] in_muL, 
            double[] in_sigmaL, double[] in_cL)
    {
        int i;
        muL = new double[L];
        sigmaL = new double[L];
        cL = new double[L];
        for(i=0; i<L; i++)
        {
            muL[i] = in_muL[i];
            sigmaL[i] = in_sigmaL[i];
            cL[i] = in_cL[i];
        }
    }
    
    public void cal_locfdr()
    {
        DRand rand_eng = new DRand();
        Normal norm_rv=new Normal(0, 1, rand_eng);

        int i,j;
        double f0_prod, f1_prod, f_prod;
        f0_prod=1; f1_prod=1; f_prod=1;
        
        for(i=0; i<group_size; i++)
        {
            f0x[i] = norm_rv.pdf(test_stat[i]);
            f1x[i] = 0;
            for(j=0;j<L; j++)
                f1x[i] = f1x[i] + cL[j] * norm_rv.pdf((test_stat[i]-muL[j])/sigmaL[j])/sigmaL[j];
            fx[i] = pi20 * f0x[i] + pi21 * f1x[i];
            
            f0_prod = f0_prod * 10 * f0x[i];
            f1_prod = f1_prod * 10 * f1x[i];
            f_prod = f_prod * 10 * fx[i];
        }
        Lfdr_g = pi0 * f0_prod / ( pi0*f0_prod +pi1* ( f_prod - Math.pow(pi20, group_size) *f0_prod )/(1-
                Math.pow(pi20, group_size ) )  );
        
        double[] num = new double[group_size];
        double den;
        den = f_prod - Math.pow(pi20, group_size) * f0_prod;
        for(i=0; i<group_size; i++ )
        {
            num[i] = pi20 * f0x[i]/fx[i] * f_prod - Math.pow(pi20, group_size)*f0_prod;
            Lfdr_gj[i] = num[i]/den;
        }
        // The following code calculate the quantities P(theta_g=1, theta_gi=1, m_gi=l|data)
        double [] each_component = new double [L];
        double total_component=0;
        double prob_theta_g_theta_gj;
        
        for( i=0; i<group_size; i++)
        {
            prob_theta_g_theta_gj = (1-Lfdr_g)*(1-Lfdr_gj[i]);
            total_component=0;
            for(j=0; j<L; j++)
            {
                each_component[j] = cL[j] * norm_rv.pdf( (test_stat[i] - muL[j])/sigmaL[j] )/sigmaL[j];
                total_component = total_component + each_component[j];
            }
            for(j=0; j<L; j++)
                prob_m_il[i][j] = prob_theta_g_theta_gj * each_component[j]/total_component;            
        }
    }
    
}
    
        