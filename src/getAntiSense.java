import java.io.*;
import java.util.*;
import java.lang.*;

public class getAntiSense
{
    public static void main(String[] args) throws IOException
    {
        String sense = "";
        String complement = "";
        String antisense = "";
        
        // Create a Scanner to read source file and convert it to an Oligo object
        Scanner inputFile = new Scanner(new File("methylation/test.txt"));
        
        //output results to a (temporary) complement file
        File fileComplement = new File("methylation/complement.txt");

        PrintWriter outputComplement = new PrintWriter(fileComplement);
        
        //obtain (temporary) complement file
        while(inputFile.hasNext())
        {
            String nt = inputFile.next();
            if(nt == "a")
                outputComplement.print("t");
            else if(nt == "c")
                outputComplement.print("g");
            else if(nt == "g")
                outputComplement.print("c");
            else if(nt == "t")
                outputComplement.print("a");
        }
        
        System.out.println("Finished creating complement file.");
        
        //no longer needed
        inputFile.close();
        outputComplement.close();
        
        // Create a Scanner to read complement file
        Scanner inputComplement = new Scanner(new File("methylation/complement.txt"));
        
        File complementFile = new File("methylation/test_AS(new).txt");
        PrintWriter fileAS = new PrintWriter(complementFile);
        
        for(long i = complementFile.length() - 1; i >= 0; i--)
            fileAS.append(inputComplement.next());
        
        System.out.println("Finished writing antisense file");
        inputComplement.close();
        fileAS.close();
    } //end main()
}
