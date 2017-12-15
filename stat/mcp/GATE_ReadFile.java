/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat.mcp;

import au.com.bytecode.opencsv.CSVReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author zhaozhg
 */
public class GATE_ReadFile {
    int G;
    int[] group_sizes;
    double[] data_all;
    String filename;
    public GATE_ReadFile(String ini_filename) throws FileNotFoundException
    {
        filename=ini_filename;
        read_file();
    }
    
    public int get_G(){ return(G); }
    public int[] get_group_sizes(){ return(group_sizes);}
    public double[] get_data_all(){ return(data_all);}
    
    private void read_file() throws FileNotFoundException
    {    
        String[] line;
        int total_hypotheses;

        try {
            CSVReader reader1 = new CSVReader(new FileReader(filename), ',');
            line = reader1.readNext();
            G = 1;
            total_hypotheses = line.length;

            while ((line = reader1.readNext()) != null) {
                G = G + 1;
                total_hypotheses = total_hypotheses + line.length;
            }
            reader1.close();

            data_all = new double[total_hypotheses];
            group_sizes = new int[G];

            int cur_ind = 0;
            int i, j, g = 0;

            reader1 = new CSVReader(new FileReader(filename), ',');

            while ((line = reader1.readNext()) != null) {
                for (i = 0; i < line.length; i++) {
                    data_all[cur_ind] = Double.parseDouble(line[i]);
                    cur_ind = cur_ind + 1;
                }
                group_sizes[g] = line.length;
                g = g + 1;
            }

        } catch (IOException | NumberFormatException e) {
            System.err.println(e);
        }
    }
    
}
