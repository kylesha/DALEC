//program will take raw Elim files and (1) reformat tag and (2) strip N's and write to a single FASTA file
import java.io.*;
import java.util.*;
import java.util.regex.*;

public class processAllFilesInDirToFASTA
{
    public static void main(String[] args) throws IOException
    {
        //create input directory object and obtain list of files in directory
        File din = new File("elim/AF104844");
        File[] list = din.listFiles();
        
        //create output file
        PrintWriter fout = new PrintWriter(new File("elim/AF104844.fa"));
        
        //read data from each of input file in the directory
        for(int i = 0; i <= list.length - 1; i++)
        {
            Scanner fin = new Scanner(list[i]);
            fin.useDelimiter("\n");
            
            //first line is id
            String ID = fin.next();
            
            //rest of file is the sequence
            String Seq = new String();
            while(fin.hasNext())
                Seq += fin.next();
            
            Pattern patID = Pattern.compile("AF.+_");
            Matcher matID = patID.matcher(ID);
            
            Pattern patSeq = Pattern.compile("N{6,}");
            Matcher matSeq = patSeq.matcher(Seq);
            
            //print id to file
            if(matID.find())
                fout.println(">" + ID.substring(matID.start(), matID.end() - 1));
            else if(! matID.find())
                fout.println(ID);
            
            //print sequence to file
            if(matSeq.find())
                fout.println(Seq.substring(0, matSeq.start()));
            else if(! matSeq.find())
                fout.println(Seq);
            
            //close current input file
            fin.close();
            
        } //end for(() loop
        fout.close(); //close output file
        System.out.println("Finished parsing!!");
    } //end main()
}
