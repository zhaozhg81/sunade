/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat;
import java.io.FileNotFoundException;
import java.io.FileReader;
import au.com.bytecode.opencsv.CSVReader;
import cern.colt.list.DoubleArrayList;
import java.io.IOException;
import sunade.stat.mcp.MCP;
import sunade.stat.mcp.sunade_MCP_Exception;

        
/**
 *
 * @author zhaozhg
 */
public class Stat {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException, sunade_MCP_Exception, sunade.stat.mcp.sunade_MCP_Exception {
        // TODO code application logic here

        int i, j, dimension=1, sample_size=1;
        double[][] dataset;
        
        DoubleArrayList Data;
        
        // Data= new DoubleArrayList(data);
        
        String datafile="data.csv";
        
        String[] line;
        try{
            CSVReader reader1 = new CSVReader( new FileReader(datafile));
            dimension=1;
            line = reader1.readNext();
            sample_size = line.length;
            
            while( reader1.readNext()!=null )
            {
                dimension=dimension+1;
            }
            reader1.close();
            
            dataset = new double[dimension][sample_size];
            
            reader1 = new CSVReader( new FileReader(datafile));
            int line_count=0;
            while( (line=reader1.readNext())!= null )
            {
                for( i=0; i< sample_size; i++)
                    dataset[line_count][i] = Double.parseDouble(line[i]);
                line_count++;
            }    
            System.out.printf("load data.");
            
            double[] data_temp;
            data_temp =new double[dimension];
            for(j=0; j<dimension;j++)
                   data_temp[j]=dataset[j][0];
            Data=new DoubleArrayList( data_temp);
            
            MCP mcp;
            String Is_P_Value;
            Is_P_Value = "test stat";
        
            mcp = new MCP(Data.copy(), Is_P_Value, 0.2);
                       
            mcp.run_analysis("locfdr");
            for(i=0;i< dimension;i++)
            {
                if( mcp.decision.get(i) )
                    System.out.printf("The %d-th hypothesis is rejected. The p-value is %.10f.\n", i, Data.get(i) );
            }

            
            
        }catch(FileNotFoundException e){
            System.out.println("File does not exist.");
        }

        
    }
    
}
