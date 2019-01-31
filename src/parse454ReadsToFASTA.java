import java.io.*;
import java.lang.*;

class parse454ReadsToFASTA
{
    public static void main(String[] args) throws IOException
    {
        File din = new File("~/DALEC/454Reads091207/");
        File[] fileList = din.listFiles();
        
        for(File file: fileList)
        {
            InputStream sin = new BufferedInputStream(new FileInputStream(file));
            PushbackInputStream fin = new PushbackInputStream(sin);
            
            String fileName = file.getName();
            
            OutputStream sout = new BufferedOutputStream(new FileOutputStream(new File("~/DALEC/454readsOUTPUT/" + fileName.substring(0, fileName.length() - 13) + ".fa")));
            PrintStream fout = new PrintStream(sout);
            
            System.out.println("Processing " + fileName + "...");
            
            while(fin.available() > 0)
            {
                char ch;
                String id = "";
                String seq = "";
                if((ch = (char)fin.read()) == '>')
                {
                    do  //extract the id tag
                    {
                        id += ch;
                    } while((ch = (char)fin.read()) != '\r' && ch != '\n');
                    
                    do  //extract the sequence readout
                    {
                        ch = (char)fin.read();
                        if(ch == 'c' || ch == 'g' || ch == 't' || ch == 'a' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'A')
                            //fout.write(ch);
                            seq += ch;
                    } while(ch != '>' && fin.available() > 0);
                    
                    fin.unread(ch);
                    fout.print(id + "\r");
                    fout.print(seq + "\r");
                } //end if
            } //end while loop
            
            fout.close();
            fin.close();
        }//end for loop
    } //end main()
}
