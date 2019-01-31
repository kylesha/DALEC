import java.io.*;

public class DpnIFrag implements Serializable
{
    private Oligo frag;     //Dpn I fragment
    private String chr;     //chromosome number
    private char strand;    //orientation(+/-)
    private int key;        //fragment number (i.e. key)
    private int start;      //start of nthFrag in genome
    private int stop;       //stop nucleotide of nthFrag
    private int fStart;     //start of fragment in fake ABI genome
    private int fStop;      //stop of fragment in fake ABI genome
    private int GATCstart;  //start of GATC in fragment
    
    //========================= constructor ==========================//
    public DpnIFrag(int key, Oligo frag, String chr, int GATCstart, int start, int stop, char strand, int fStart, int fStop)
    {
        this.key = key;
        this.frag = frag;
        this.chr = chr;
        this.GATCstart = GATCstart;
        this.start = start;
        this.stop = stop;
        this.strand = strand;
        this.fStart = fStart;
        this.fStop = fStop;
    }
    
    //empty constructor
    public DpnIFrag()
    {
//        this.key = 0;
//        this.frag = new String();
//        this.chr = "";
//        this.GATCstart = 0;
//        this.start = 0;
//        this.stop = 0;
//        this.strand = ' ';
//        this.fStart = 0;
//        this.fStop = 0;
    }
    
    
    //====================== methods ==========================//
    public String toString()
    {
        return (this.key - 9) + '\t' + this.frag.toString() + '\t' + this.chr + '\t' + this.GATCstart + '\t' + this.start + '\t' + this.stop + '\t' + this.strand + '\t' + this.fStart + '\t' + this.fStop;
    }
    
    
    public String getChr()
    {
        return this.chr;
    }
    
    
    public int getfStart()
    {
        return this.fStart;
    }
    
    
    public int getfStop()
    {
        return this.fStop;
    }
    
    
    public int getGATCstart()
    {
        return GATCstart;
    }
    
    
    public int getKey()
    {
        return this.key;
    }
    
    
    public Oligo getSeq()
    {
        return this.frag;
    }
    
    
    public int getStart()
    {
        return this.start;
    }
    
    
    public int getStop()
    {
        return this.stop;
    }
    
    
    public char getStrand()
    {
        return this.strand;
    }
    
    
    public static DpnIFrag toDpnIObj(String str)
    {
        String[] input = str.split("\t");
        
        int key = Integer.parseInt(input[0]);
        Oligo frag = new Oligo(input[1]);
        String chr = input[2];
        int GATCstart = Integer.parseInt(input[3]);
        int start = Integer.parseInt(input[4]);
        int stop = Integer.parseInt(input[5]);
        char strand = input[6].charAt(0);
        int fStart = Integer.parseInt(input[7]);
        int fStop = Integer.parseInt(input[8]);
        
        return new DpnIFrag(key, frag, chr, GATCstart, start, stop, strand, fStart, fStop);
    }
}

