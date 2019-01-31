import java.io.*;
import java.util.*;
import java.lang.*;

public class InSilicoCloning
{
    public static void main(String[] args) throws Exception
    {
        Oligo target = new Oligo("gatc");
        String sequence = "";
        int allowedMismatches = 0;
        
        // Create a Scanner to read source file and convert it to an Oligo object
        Scanner inputFile = new Scanner(new File("methylation/WS140Chrom1.txt"));
        
        //convert file object into a sequence
        while(inputFile.hasNext())
        {
            sequence += inputFile.next();
        }
        
        //convert sequence to an Oligo
        Oligo chrom1 = new Oligo(sequence);
        ArrayList list = chrom1.getAllMatchCoordinates(target, allowedMismatches);
        
        //output results to a txt file
        File outputFile = new File("methylation/WS140Chrom1_output.txt");
        PrintWriter output = new PrintWriter(outputFile);
        
        output.println("LG\tPOSITION\tSENSE\t\t\tANTISENSE");
        for(int i = 0; i <= list.size() - 1; i++)
        {
            int start = Integer.parseInt((list.get(i)).toString()) + 2;
            int end = start + 20;
            Oligo tempOligo = new Oligo(chrom1.extractSequence(start, end)); //extract sequence
            Oligo tempOligoAp = new Oligo(tempOligo.antiparallel()); //extract antiparallel sequence

            //write to file
            output.println("I\t" + start + "\t\t" + tempOligo.toString() + "\t" + tempOligoAp.toString());
        }
        
        // Close the files
        inputFile.close();
        output.close();
        System.out.println("Done!!! Size of list is: " + list.size());
    } //end main()
}


