import java.io.*;
import java.util.regex.*;

class WS170InSilicoClonesABIv3
{
    public static void main(String[] args) throws IOException
    {
        //create input directory object and obtain list of files in directory
        File din = new File("~/DALEC/WS170_input/");
        File[] list = din.listFiles();
        
        if(din.exists())
        {
            OutputStream real = new BufferedOutputStream(new FileOutputStream("~/DALEC/WS170_in_silico_clones/WS170_InSilicoClones.txt"));
            PrintStream ABIreal = new PrintStream(real);
            OutputStream fake = new BufferedOutputStream(new FileOutputStream("~/DALEC/WS170_in_silico_clones/WS170_fake_genome.txt"));
            PrintStream ABIfake = new PrintStream(fake);
            
            int nthFrag = 0;            //Frag number
            int fakeStart = 0;          //start coordinate of nth Frag in fake ABIreal
            int fakeStop = 0;           //stop coordinate of nth Frag in fake ABIreal
            final int bufSize = 60;
            String linker = new String("gtcggactgg");   //10 nucleotide constant linker region
            
            Pattern pat = Pattern.compile("[acgt]{28}gatc[acgt]{28}");    //regex pattern to search for
            
            //cycle through all files in directory
            for(File inputFile: list)
            {
                InputStream in = new BufferedInputStream(new FileInputStream(inputFile));
                PushbackInputStream fin = new PushbackInputStream(in, bufSize);
                
                char ch;
                int senseStart = 0;                 //start coordinate of sense fragment
                int senseStop = 0;                  //stop coordinate of sense fragment
                int senseGATC = 0;                  //start of 'g' in "gatc" of strand fragment
                int antisenseStart = 0;             //start coordinate of antisense fragment
                int antisenseStop = 0;              //stop coordinate of antisense fragment
                int antisenseGATC = 0;              //start of 'g' in "gatc" of antisense fragment
                byte[] buffer = new byte[bufSize];  //sliding window
                String strand = "";                 //orientation (+/-) strand
                String id = "";                     //FASTA id
                String chr = "";                    //chromosome number
                String window = "";                 //sliding window
                String senseFrag = "";              //found DpnI Frag on sense strand
                String antisenseFrag = "";          //found DpnI Frag on antisense strand
                
                if((ch = (char)fin.read()) == '>')
                {
                    //extract the chromosome number
                    while((ch = (char)fin.read()) != '\r')
                        id += ch;
                    chr = id.substring(0);
                    
                    System.out.println("FASTA id found. Processing " + inputFile + "...");
                    System.out.println();
                    
                    //find (+/-)gatc frags in current chromosome
                    for(int i = 0; fin.available() >= bufSize; i++)
                    {
                        fin.read(buffer);
                        
                        if(buffer[28] == 'g')
                        {
                            window = new String(buffer);
                            Matcher mat = pat.matcher(window.toLowerCase());
                            
                            if(mat.matches())
                            {
                                //process sense fragment
                                strand = "+";
                                nthFrag++;
                                senseStart = i;
                                senseStop = i + 29;
                                senseGATC = senseStop - 1;                //start of 'g' in "[n]28gatc" on sense strand
                                fakeStart = 40 * (nthFrag - 1) + 1; //start coordinate of nth Frag in fake genome
                                fakeStop = 40 * nthFrag;            //stop coordinate of nth Frag in fake genome
                                senseFrag = window.substring(0, 30);
                                ABIreal.print(nthFrag + "\t" + senseFrag + "\t" + chr + "\t" + senseGATC + "\t" + senseStart + "\t" + senseStop + "\t" + strand + "\t" + fakeStart + "\t" + fakeStop + "\r");
                                ABIfake.print(senseFrag + linker + "\n");
                                
                                //process antisense fragment
                                strand = "-";
                                nthFrag++;
                                antisenseStart = senseStop + 1;
                                antisenseStop = senseStart + 59;
                                antisenseGATC = senseGATC + 3;
                                fakeStart = 40 * (nthFrag - 1) + 1; //start coordinate of nth Frag in fake genome
                                fakeStop = 40 * nthFrag;            //stop coordinate of nth Frag in fake genome
                                antisenseFrag = new Oligo(window.substring(30)).antiparallel().toString();  //31st nucleotide = 30th character of array
                                ABIreal.print(nthFrag + "\t" + antisenseFrag + "\t" + chr + "\t" + antisenseGATC + "\t" + antisenseStart + "\t" + antisenseStop + "\t" + strand + "\t" + fakeStart + "\t" + fakeStop + "\r");
                                ABIfake.print(antisenseFrag + linker + "\n");
                            }
                        }
                        
                        fin.unread(buffer, 1, bufSize - 1); //pushback all characters in buffer except the first
                    } //end for loop
                } //end if((ch = (char)fin.read()) == '>')
                else System.out.println("File tag not found. Skipping file..." + inputFile);
                
                fin.close();
            } //end outer for loop
            
            ABIreal.close();
            ABIfake.close();
        } //end if(din.exists())
        else
        {
            System.out.println(din + "does not exist. Program terminated");
            System.exit(0);
        }
        
        System.out.println("Finished processing entire directory.");
    } //end main()
}






