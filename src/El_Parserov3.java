/**
 *@author Ky Sha
 *@version 1.7(07/11/07) - parse5PrimeLinker and parse3PrimeLinker revised to incorporate fuzzySmithWaterman() into the search strategy
 *@version 1.8(07/16/07) - output modified to tabulated format
 *@version 1.9(07/23/07) - parses 454, ABI, and Solexa readouts and outputs to each respective file
 *@version 2.0(08/01/07) - incorporates exception handling
 *@version 3.0(08/17/07) - completely re-written algorithm to use the more intelligent combined fuzzyMatch/smithWaterman methods
 *@version 3.1(08/21/07) - parse5PrimeLinker() modified to more intelligently locate MmeI site if a perfect match doesn't exist
 *@version 3.2(08/23/07) - renaming of some variables to make them unambiguous. Oligos shown in 5'-->3' orientation.
 */
import java.io.*;
import java.lang.*;
import java.util.ArrayList;
import java.util.regex.*;
import java.util.Scanner;

class El_Parserov3_2
{
    public static void main(String[] args) throws IOException
    {
        File file = new File("/Users/kylesha/Desktop/test.txt");
        PrintWriter fout = new PrintWriter(new File("/Users/kylesha/Desktop/test_parsed01.txt"));
        Scanner fin = new Scanner(file);
        
        String id = "";             //id of oligo
        String orientation = "";    //orientation of oligo (SENSE/ANTISENSE)
        String exp = "";            //experiment number
        String line = "";           //worm line used
        String comments = "";       //AB-, A-, or B-linker successfully parsed
        String platform = "";       //sequencing platform
        Oligo oligo = new Oligo();
        Oligo parsed5Prime = new Oligo();
        Oligo parsed3Prime = new Oligo();
        Oligo parsedBoth = new Oligo();
        Oligo insert = new Oligo();         //both_parsed/left_parsed/right_parsed/unparsed oligo
        Oligo _5primeLinker = new Oligo();
        Oligo _3primeLinker = new Oligo();
        Oligo _5primeQuery = new Oligo();
        Oligo _3primeQuery = new Oligo();
        Oligo _5primeLinker454A1_1 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATTCTCTCCGAC");            //454 5' linker
        Oligo _3primeLinker454B1_3 = new Oligo("ATCACCAGCTGCTCATAGACTGTACCAGTCTGAGCGGGCTGGCAAGGC");            //454 3' linker
        Oligo _5primeLinkerABIf1a = new Oligo("AACTGCCCCGGGTTCCTCATTCTCTGGACCTGCTCGATCCAGTCCGAC");             //ABI 5' linker
        Oligo _3primeLinkerABIf2a = new Oligo("ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGGTTGACTA");             //ABI 3' linker
        Oligo _5primeLinkerSolexA1_1 = new Oligo("CAAGCAGAAGACGGCATACGATCCTGAGTACACTATGTTCCGAC");             //Solexa 5' linker
        Oligo _3primeLinkerSolexB1_1P = new Oligo("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCGGTGGTCGCCGTATCATT");    //Solexa 3' linker
        Oligo _5primeQuery454A = new Oligo("GTTCAAGTCGTTCCTCATTCTC");       //454 5' query
        Oligo _3primeQuery454B = new Oligo("ATCACCAGCTGCTCATAGA");          //454 3' query
        Oligo _5primeQueryABI1 = new Oligo("CATTCTCTGGACCTGCTCGATCCAG");    //ABI 5' query
        Oligo _3primeQueryABI2 = new Oligo("ATCACCGACTGCCCATAGAG");         //ABI 3' query
        Oligo _5primeQuerySolexA1_1 = new Oligo("CGATCCTGAGTACACTATGT");    //Solexa 5' query
        Oligo _3primeQuerySolexB1_1P = new Oligo("AGATCGGAAGAGCGTCGTGTA");  //Solexa 3' query
        Oligo MmeISite = new Oligo("TCCGAC");
        int min5primeKeyLength = _5primeQuery.length() - 2;;
        int min3primeKeyLength = _3primeQuery.length() - 2;;
        final int MIN_LENGTH = 200; //minimum length for readout to be processed
        int allowedMismatches = 3;
        int ins = 2;
        int del = 2;
        
//        Oligo A02 = new Oligo   ("TGTTGGAGGGGTTCGTCTT-ATGATACGGCGACCACCGACACTCTTTCCCTACACGACGCTCTTCCGATCTATAATTTTCCTCACCCTGAGTCGGAACATAGTGTACTCAGGATCGTATGCCGTCTTCTGCTTGAAGGGCGAATTCGTTTAAACCTGCAGGACTAGTCCCTTTAGTGAGGGTTAATTCTGAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGA");
//        Oligo B12 = new Oligo("AGCACGCTCCGGATTCGCCCTTAATGATACGGCGACCACCGACACTCTTTCCCTACACGACGCTCTTCCGATCTATTCGCCTCCTGTGATGTTGAGTCGGAACATAGTGTACTCAGGATCGTATGCCGTCTTCTGCTTGAAGGGCGAATTCGTTTAAACCTGCAGGACTAGTCCCTTTAGTGAGGGTTAATTCTGAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGA");
//        Oligo F10 = new Oligo("TTCTTAGTCCGGATTCGCCTTCAGCAGAAGACGGCATACGATCCTGAGTACACTATGTTCCGACCTATGTTCCGACGTCGGAACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCGGTGGTCGCCGTATCATTAAGGGCGAATTCGTTTAAACCTGCAGGACTAGTCCCTTTAGTGAGGGTTAATTCTGAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGA");
//
//        try
//        {
//            System.out.println("found5PrimeLinker(F10): " + found5PrimeLinker(F10, _5primeLinkerSolexA1_1));
//            System.out.println("found3PrimeLinker(F10): " + found3PrimeLinker(F10,_3primeLinkerSolexB1_1P));
//            System.out.println("getSense(F10): " + getSense(F10, _5primeLinkerSolexA1_1, _3primeLinkerSolexB1_1P));
//        }
//        catch (OligoException error)
//        {System.out.println(error);}
        
        while(fin.hasNext())
        {
            try
            {
                id = fin.next();  //read id
                oligo = new Oligo(fin.next()); //read corresponding oligo
                
                //determine which platform and which linkers to parse
                //begin 454
                if(isPlatform(oligo, _5primeLinker454A1_1, _3primeLinker454B1_3))
                {
                    platform = "454";
                    exp = "11th";
                    orientation = getSense(oligo, _5primeLinker454A1_1, _3primeLinker454B1_3);
                    
                    if(orientation == "SENSE")
                    {
                        _5primeLinker = _5primeLinker454A1_1;
                        _3primeLinker = _3primeLinker454B1_3;
                        _5primeQuery = _5primeQuery454A;
                        _3primeQuery = _3primeQuery454B;
                    }
                    else if(orientation == "SENSE_5prime")
                    {
                        _5primeLinker = _5primeLinker454A1_1;
                        _5primeQuery = _5primeQuery454A;
                    }
                    else if(orientation == "SENSE_3prime")
                    {
                        _3primeLinker = _3primeLinker454B1_3;
                        _3primeQuery = _3primeQuery454B;
                    }
                    else if(orientation == "ANTISENSE")
                    {
                        _5primeLinker = _3primeLinker454B1_3.antiparallel();
                        _3primeLinker = _5primeLinker454A1_1.antiparallel();
                        _5primeQuery = _3primeQuery454B.antiparallel();
                        _3primeQuery = _5primeQuery454A.antiparallel();
                    }
                    else if(orientation == "ANTISENSE_5prime")
                    {
                        _5primeLinker = _3primeLinker454B1_3.antiparallel();
                        _5primeQuery = _3primeQuery454B.antiparallel();
                    }
                    else if(orientation == "ANTISENSE_3prime")
                    {
                        _3primeLinker = _5primeLinker454A1_1.antiparallel();
                        _3primeQuery = _5primeQuery454A.antiparallel();
                    }
                } //end 454
                
                //begin ABI
                else if(isPlatform(oligo, _5primeLinkerABIf1a, _3primeLinkerABIf2a))
                {
                    System.out.println("ABI");
                    platform = "ABI";
                    exp = "15th";
                    orientation = getSense(oligo, _5primeLinkerABIf1a, _3primeLinkerABIf2a);
                    
                    if(orientation == "SENSE")
                    {
                        _5primeLinker = _5primeLinkerABIf1a;
                        _3primeLinker = _3primeLinkerABIf2a;
                        _5primeQuery = _5primeQueryABI1;
                        _3primeQuery = _3primeQueryABI2;
                    }
                    else if(orientation == "SENSE_5prime")
                    {
                        _5primeLinker = _5primeLinkerABIf1a;
                        _5primeQuery = _5primeQueryABI1;
                    }
                    else if(orientation == "SENSE_3prime")
                    {
                        _3primeLinker = _3primeLinkerABIf2a;
                        _3primeQuery = _3primeQueryABI2;
                    }
                    else if(orientation == "ANTISENSE")
                    {
                        _5primeLinker = _3primeLinkerABIf2a.antiparallel();
                        _3primeLinker = _5primeLinkerABIf1a.antiparallel();
                        _5primeQuery = _3primeQueryABI2.antiparallel();
                        _3primeQuery = _5primeQueryABI1.antiparallel();
                    }
                    else if(orientation == "ANTISENSE_5prime")
                    {
                        _5primeLinker = _3primeLinkerABIf2a.antiparallel();
                        _5primeQuery = _3primeQueryABI2.antiparallel();
                    }
                    else if(orientation == "ANTISENSE_3prime")
                    {
                        _3primeLinker = _5primeLinkerABIf1a.antiparallel();
                        _3primeQuery = _5primeQueryABI1.antiparallel();
                    }
                } //end ABI
                
                //begin Solexa
                else if(isPlatform(oligo, _5primeLinkerSolexA1_1, _3primeLinkerSolexB1_1P))
                {
                    System.out.println("SOLEXA");
                    platform = "SOLEXA";
                    exp = "16th";
                    orientation = getSense(oligo, _5primeLinkerSolexA1_1, _3primeLinkerSolexB1_1P);
                    
                    if(orientation == "SENSE")
                    {
                        _5primeLinker = _5primeLinkerSolexA1_1;
                        _3primeLinker = _3primeLinkerSolexB1_1P;
                        _5primeQuery = _5primeQuerySolexA1_1;
                        _3primeQuery = _3primeQuerySolexB1_1P;
                    }
                    else if(orientation == "SENSE_5prime")
                    {
                        _5primeLinker = _5primeLinkerSolexA1_1;
                        _5primeQuery = _5primeQuerySolexA1_1;
                    }
                    else if(orientation == "SENSE_3prime")
                    {
                        _3primeLinker = _3primeLinkerSolexB1_1P;
                        _3primeQuery = _3primeQuerySolexB1_1P;
                    }
                    else if(orientation == "ANTISENSE")
                    {
                        _5primeLinker = _3primeLinkerSolexB1_1P.antiparallel();
                        _3primeLinker = _5primeLinkerSolexA1_1.antiparallel();
                        _5primeQuery = _3primeQuerySolexB1_1P.antiparallel();
                        _3primeQuery = _5primeQuerySolexA1_1.antiparallel();
                    }
                    else if(orientation == "ANTISENSE_5prime")
                    {
                        _5primeLinker = _3primeLinkerSolexB1_1P.antiparallel();
                        _5primeQuery = _3primeQuerySolexB1_1P.antiparallel();
                    }
                    else if(orientation == "ANTISENSE_3prime")
                    {
                        _3primeLinker = _5primeLinkerSolexA1_1.antiparallel();
                        _3primeQuery = _5primeQuerySolexA1_1.antiparallel();
                    }
                }//end Solexa
                
                //actual parsing process
                if(oligo.isDNA() && oligo.length() >= MIN_LENGTH && oligo.countChar('n') < 5)
                {
                    if(orientation == "SENSE")
                    {
                        parsed5Prime = parse5PrimeLinker(oligo, _5primeLinker, _5primeQuery, MmeISite, allowedMismatches, ins, del, min5primeKeyLength);
                        parsedBoth = parse3PrimeLinker(parsed5Prime, _3primeLinker, _3primeQuery, allowedMismatches, ins, del, min3primeKeyLength);
                        insert = parsedBoth;
                        comments = "SENSE both parsed";
                    }
                    else if(orientation == "ANTISENSE")
                    {
                        parsed5Prime = parse5PrimeLinker(oligo, _5primeLinker, _5primeQuery, MmeISite.antiparallel(), allowedMismatches, ins, del, min5primeKeyLength);
                        parsedBoth = parse3PrimeLinker(parsed5Prime, _3primeLinker, _3primeQuery, allowedMismatches, ins, del, min3primeKeyLength);
                        insert = parsedBoth.reverse();
                        comments = "ANTISENSE both parsed";
                    }
                    else if(orientation == "SENSE_5prime")
                    {
                        parsed5Prime = parse5PrimeLinker(oligo, _5primeLinker, _5primeQuery, MmeISite, allowedMismatches, ins, del, min5primeKeyLength);
                        insert = parsed5Prime;
                        comments = "SENSE_5prime parsed";
                    }
                    else if(orientation == "SENSE_3prime")
                    {
                        parsed3Prime = parse3PrimeLinker(oligo, _3primeLinker, _3primeQuery, allowedMismatches, ins, del, min3primeKeyLength);
                        insert = parsed3Prime;
                        comments = "SENSE_3prime parsed";
                    }
                    else if(orientation == "ANTISENSE_5prime")
                    {
                        parsed5Prime = parse5PrimeLinker(oligo, _5primeLinker, _5primeQuery, MmeISite, allowedMismatches, ins, del, min5primeKeyLength);
                        insert = parsed5Prime;
                        comments = "ANTISENSE_5prime parsed";
                    }
                    else if(orientation == "ANTISENSE_3prime")
                    {
                        parsed3Prime = parse3PrimeLinker(oligo, _3primeLinker, _3primeQuery, allowedMismatches, ins, del, min3primeKeyLength);
                        insert = parsed3Prime.reverse();
                        comments = "ANTISENSE_3prime parsed";
                    }
                    else if(orientation == "UNKNOWN")
                    {
                        insert = oligo;
                        comments = "UNKNOWN orientation";
                    }
                } //outer if()
                else
                {
                    insert = oligo;
                    orientation = "N/A";
                    comments = "FAILED: insert to short or too many N's";
                    fout.println(id + '\t' + insert + '\t' + orientation + '\t' + platform + '\t' + exp + '\t' + comments);
                }
                fout.println(id + '\t' + insert + '\t' + orientation + '\t' + platform + '\t' + exp + '\t' + comments);
            } //end try{}
            catch (OligoException error)
            {
                insert = oligo;
                orientation = "N/A";
                fout.println(id + '\t' + insert + '\t' + orientation + '\t' + platform + '\t' + exp + '\t' + error);
            }
        } //end while() loop
        fin.close();
        fout.close();
    }//end main()
    
    
    
    //============================= El Parsero accessory methods ==============================//
    /**
     *Assumes insert is in 5'-->3' orientation
     *Parses out the left linker if anchor insert (MmeIsite) is found. Will parse out entire left linker including the MmeI site
     */
    private static Oligo parse5PrimeLinker(Oligo oligo, Oligo _5primeLinker, Oligo query, Oligo anchor, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        int startIndex = _5primeLinker.toString().indexOf(query.toString());
        int endIndex = startIndex + query.length() - 1;
        int anchorIndex = _5primeLinker.toString().indexOf(anchor.toString());
        int allowedAnchorMismatch = 0;
        final int LEFTLINKER_LENGTH = _5primeLinker.length();
        boolean found = false;
        boolean anchorExists = false;
        Oligo anchor2 = new Oligo(_5primeLinker.toString().substring(anchorIndex - 1));
        Oligo slidingKey = new Oligo();
        Oligo foundKey = new Oligo();
        Oligo swKey = new Oligo();
        
        //if perfect match not found for MmeI site, extend anchor to include 1nt left and allow 1 mismatch
        if(oligo.isFuzzyMatch(anchor, 0))
            anchorExists = true;
        else if(oligo.isFuzzyMatch(anchor2, 1))
        {
            anchorExists = true;
            anchor = anchor2;
            allowedAnchorMismatch = 1;
        }
        else anchorExists = false;
        
        if(anchorExists)
        {
            final int ANCHOR_START = oligo.getFirstMatchCoordinate(anchor, allowedAnchorMismatch);
            final int ANCHOR_STOP = ANCHOR_START + anchor.length() - 1;
            Oligo searchWindow = new Oligo(oligo.extractSequence(0, ANCHOR_START - 1)); //limit search window to left of MmeI site
            
            for(int i = 0; startIndex > i && !found; i++) //shift slidingKey one nt left until reach 5' of left linker
            {
                try
                {
                    slidingKey = _5primeLinker.extractSequence(startIndex - i, endIndex - i);
                    swKey = searchWindow.getFuzzySWkey(slidingKey, mismatches, ins, del, minKeyLength);
                    foundKey = oligo.extractSequence(oligo.getFirstMatchCoordinate(swKey, 0), ANCHOR_STOP);
                    if(oligo.isFuzzyMatch(foundKey, 0))
                        found = true;
                }
                catch (OligoException error)
                {continue;}
            }
            
            if(found)
                return oligo.exciseLeftFrom(foundKey, 0);
            else throw new OligoException("MatchNotFound ", "parse5PrimeLinker()");
        }
        else throw new OligoException("MmeI site not found", "parse5PrimeLinker(): ");
    } //end parse5PrimeLinker() method
    
    
    
    private static Oligo parse3PrimeLinker(Oligo oligo, Oligo _3primeLinker, Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        int startIndex = _3primeLinker.toString().indexOf(query.toString());
        int endIndex = startIndex + query.length() - 1;
        final int RIGHTLINKER_LENGTH = _3primeLinker.length();
        Oligo swKey = new Oligo();
        Oligo foundKey = new Oligo();
        
        try
        {
            if(oligo.isFuzzySmithWatermanMatch(query, mismatches, ins, del, minKeyLength))
            {
                swKey = oligo.getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
                foundKey = oligo.extractSequence(swKey, mismatches);
                return oligo.exciseRightFrom(foundKey,0);
            }
            else throw new OligoException("MatchNotFound", "parse3PrimeLinker()");
        }
        catch (OligoException error)
        {throw error;}
    } //end parse3PrimeLinker()
    
    
    
    //---------------------------------- getSense() -----------------------------------//
    private static String getSense(Oligo oligo, Oligo _5primeLinker, Oligo _3primeLinker) throws OligoException
    {
        String sense = new String();
        
        try
        {
            boolean SENSE_foundLinkerL = found5PrimeLinker(oligo, _5primeLinker);
            boolean SENSE_foundLinkerR = found3PrimeLinker(oligo, _3primeLinker);
            boolean ANTISENSE_foundLinkerL = found5PrimeLinker(oligo, _3primeLinker.antiparallel());
            boolean ANTISENSE_foundLinkerR = found3PrimeLinker(oligo, _5primeLinker.antiparallel());
            
            if(SENSE_foundLinkerL && SENSE_foundLinkerR)
                sense = "SENSE";
            else if(ANTISENSE_foundLinkerL && ANTISENSE_foundLinkerR)
                sense = "ANTISENSE";
            else if (SENSE_foundLinkerL ^ SENSE_foundLinkerR)
                sense = SENSE_foundLinkerL ? "SENSE_5prime" : "SENSE_3prime";
            else if (ANTISENSE_foundLinkerL ^ ANTISENSE_foundLinkerR)
                sense = ANTISENSE_foundLinkerL ? "ANTISENSE_5prime" : "ANTISENSE_3prime";
            else sense = "UNKNOWN";
        }
        catch (OligoException error)
        {throw error;}
        
        return sense;
    } //end getSense() method
    
    
    
    //---------------------------------- getSense() -----------------------------------//
    private static boolean isPlatform(Oligo oligo, Oligo _5primeLinker, Oligo _3primeLinker) throws OligoException
    {
        boolean platform = false;
        
        try
        {
            boolean SENSE_found5PrimeLinker = found5PrimeLinker(oligo, _5primeLinker);
            boolean SENSE_found3PrimeLinker = found3PrimeLinker(oligo, _3primeLinker);
            boolean ANTISENSE_found5PrimeLinker = found5PrimeLinker(oligo, _3primeLinker.antiparallel());
            boolean ANTISENSE_found3PrimeLinker = found3PrimeLinker(oligo, _5primeLinker.antiparallel());
            
            if((SENSE_found5PrimeLinker && SENSE_found3PrimeLinker) || (ANTISENSE_found5PrimeLinker && ANTISENSE_found3PrimeLinker))
                platform = true;
            else throw new OligoException("UNKNOWN PLATFORM");
        }
        catch (OligoException error)
        {throw error;}
        
        return platform;
    } //end getSense() method
    
    
    
    //------------------------ found5PrimeLinker() ----------------------------//
    private static boolean found5PrimeLinker(Oligo oligo, Oligo _5primeLinker) throws OligoException
    {
        try
        {
            Oligo slidingQuery = new Oligo();
            Oligo searchWindow = oligo.extractSequence(0, (oligo.length() - 1) / 2); //search only first half of the oligo
            int windowWidth = 11;
            int found = 0;
            
            for(int i = 0; (i + windowWidth < _5primeLinker.length() - 1)  &&  found <= 3; i++)
            {
                slidingQuery = _5primeLinker.extractSequence(i, i + windowWidth); //scan along oligo to get key of 11 nucleotides
                if(searchWindow.isFuzzyMatch(slidingQuery, 0))
                    ++found;
            }
            return (found >= 3) ? true : false;
        }
        catch (OligoException error)
        {throw error;}
    } //end found5PrimeLinker() method
    
    
    
    //------------------------ found3PrimeLinker() ----------------------------//
    private static boolean found3PrimeLinker(Oligo oligo, Oligo _3primeLinker) throws OligoException
    {
        try
        {
            Oligo slidingQuery = new Oligo();
            Oligo searchWindow = oligo.extractSequence(oligo.length() / 2, oligo.length() - 1); //search only second half of the oligo
            int found = 0;
            
            for(int i = _3primeLinker.length() - 1; i > 11 && found <= 3; i--)
            {
                slidingQuery = _3primeLinker.extractSequence(i - 11, i); //get key of 11 nucleotides and scan along oligo
                if(searchWindow.isFuzzyMatch(slidingQuery, 0))
                    ++found;
            }
            return (found >= 3) ? true : false;
        }
        catch (OligoException error)
        {throw error;}
    } //end found5PrimeLinker() method
}//end El Parsero

