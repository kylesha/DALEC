//takes entire directory and converts to a single FASTA file, with trailing N's stripped
import java.io.*;
import java.util.*;

public class toFASTA
{
    public static void main(String[] args) throws IOException
    {
        File din = new File("~/DALEC/454Reads091207");
        File[] list = din.listFiles();
        
        String temp = new String();
        String id = new String();
        String seq = new String();
        
        for(File file: list)
        {
            Scanner fin = new Scanner(file);
            String fileName = new String(file.getName());
            fin.useDelimiter("\r");
            
            PrintWriter fout = new PrintWriter(new File("~/DALEC/454readsOUTPUT/" + fileName.substring(0, fileName.length() - 13) + ".fa"));
            
            while(fin.hasNextLine())
            {
                for(int i = 1; i <= 3;)
                {
                    temp = fin.nextLine();
                    
                    if(temp.contains(">"))
                    {
                        id = temp;
                        i++;
                    }
                    else
                    {
                        seq += temp;
                        i++;
                    }
                }
                //write id and sequence to file
                fout.println(id);
                fout.println(seq);
                temp = "";
                id = "";
                seq = "";
            }//end while loop
            
            fin.close();
            fout.close();
        }
    } //end main()
}

