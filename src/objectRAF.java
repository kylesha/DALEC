import java.io.*;
import java.util.*;

public class objectRAF
{
    public static void main(String[] args) throws IOException
    {
        try
        {
            //open inputstream and define DpnIFrag object
            Scanner fin = new Scanner(new File("~/DALEC/test/test_input.txt"));
            DpnIFrag dpn = new DpnIFrag();
            
            //obtain byte array for temp storage of objects
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream oos = new ObjectOutputStream(bos);
            
            ArrayList objectList = new ArrayList<int[]>();
            
            RandomAccessFile raf = new RandomAccessFile("~/DALEC/test/testRAF.bin", "rw");
            
            int bufLength;
            int position = 0;
            //serialize objects to file, one at a time
            while(fin.hasNextLine())
            {
                dpn = DpnIFrag.toDpnIObj(fin.nextLine());
                oos.writeObject(dpn);
//            System.out.println("size of bos: " + bos.size());
                
                //convert byte array to bytes
                byte[] buf = bos.toByteArray();
                int[] objectInfo = new int[2];
                
                bufLength = buf.length;
                objectInfo[0] = position;    //get the identifier of present object
                objectInfo[1] = bufLength;      //get the length of present object
                objectList.add(objectInfo);
                System.out.println("position: " + objectInfo[0] + " bufLength: " + objectInfo[1]);
                System.out.println(dpn);
                
                //write to file
                raf.seek(raf.length());
                raf.write(buf);
                position += bufLength;
                bos.reset();
                oos.flush();
            }
            
            //close input stream
            fin.close();
            oos.close();
            bos.close();
            System.out.println("finished serializing");
            
            
            int index = 2;
            int[] objectInfo = (int[])objectList.get(index);
            byte[] buf = new byte[objectInfo[1]];
            
            System.out.println("objectInfo[0] = " + objectInfo[0]);
            System.out.println("objectInfo[1] = " + objectInfo[1]);
            
            raf.seek(objectInfo[0]);
            raf.readFully(buf);
            
            ByteArrayInputStream bis = new ByteArrayInputStream(buf);
            ObjectInputStream ois = new ObjectInputStream(bis);

            dpn = (DpnIFrag)ois.readObject();
            System.out.println(dpn);
            ois.close();
            
        }
        catch (FileNotFoundException ex)
        {
            ex.printStackTrace();
        }
        catch (IOException ex)
        {
            ex.printStackTrace();
        }
        catch (ClassNotFoundException ex)
        {
            ex.printStackTrace();
        }
    }//end main()
}
