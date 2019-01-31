//parse DpnI-MmeI Cloning Method in silico clones using streams instead of the Oligo class; does entire directory
import java.io.*;
import java.lang.*;
import java.util.regex.*;

class InSilicoCloningStreamsProcessAll
{
    public static void main(String[] args) throws IOException
    {
        //create input directory object and obtain list of files in directory
        File din = new File("/Users/kylesha/Desktop/test");
        File[] list = din.listFiles();
        
        for(int f = 0; f <= list.length - 1; f++) //cycle through all files in directory
        {
            InputStream fin = new BufferedInputStream(new FileInputStream(list[f]));
            String outputFile = new String("/Users/kylesha/Desktop/test_output.txt");
            PrintWriter fout = new PrintWriter(new File(outputFile));
            
            //extract the id only
            String id = new String();
            char ch = (char)fin.read();
            if(ch == '>')
            {
                do
                {
                    id += ch;
                } while((ch = (char)fin.read()) != '\r');
            }
            else
            {
                fin.close();
                System.out.println("Title tag not found. Program terminated");
                return;
            }
            
            //if last ch read is a CR, then title is found ==> proceed to parse in silico clones
            if((char)fin.read() != '\r')
            {
                System.out.println("File tag found. Processing " + list[f] + "...");
                
                String tag = new String(id.substring(1, id.length())); //tag = chromosome number and genome version
                String target = new String("gatc................");
                int n = target.length() + 2; //size of sliding window
                
                //pattern to be matched is target plus all 16 combos of any two nucleotides
                Pattern pat = Pattern.compile(target + "(aa|ac|ag|at|ca|cc|cg|ct|ga|gc|gg|gt|ta|tc|tg|tt)");
                for(int i = 1; fin.available() - 1 > n; i++) //last character is CR, don't read!!
                {
                    fin.mark(n); //mark current nucleotide index at i
                    byte[] buffer = new byte[n]; //create buffer of size n to hold sliding window
                    fin.read(buffer);
                    String window = new String(buffer);  //obtain sliding window of size n, starting at i
                    
                    //if found target match within current sliding window, send to output file
                    Matcher mat = pat.matcher(window.toLowerCase());
                    if(mat.matches())
                    {
                        fout.print(tag + "\t" + i + "\t" + mat.group());
                        fout.println();
                    }
                    
                    fin.reset(); //go back to i
                    fin.read(); //discard ith ch
                } //for loop
            }
            else //abnormal termination because tag was not found
            {
                fin.close();
                fout.close();
                System.out.println("Title tag not found! Program terminated");
                return;
            }
            
            //normal termination
            fin.close();
            fout.close();
            System.out.println("Finished parsing " + list[f] + ". Output file written to: " + outputFile);
        } //end first for loop
        
        System.out.println("Finished processing entire directory.");
    } //end main()
}

