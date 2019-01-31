/**
 *@author Ky Sha
 *@version 2.0(10/6/07)
 */
import java.io.*;
import java.util.Scanner;

class parse454v2
{
    public static void main(String[] args) throws IOException, OligoException
    {
        File din = new File("~/DALEC/454reads_input");
        File[] fileList = din.listFiles();
        
//        PrintWriter parsed = new PrintWriter(new File("~/DALEC/454reads_output/parsed_454reads.txt"));
//        PrintWriter unknown = new PrintWriter(new File("~/DALEC/454reads_output/unknown_454reads.txt"));
//        PrintWriter not_mine = new PrintWriter(new File("~/DALEC/454reads_output/not_mine_454reads.txt"));
        
        int notMine = 0;
        int unknown_seq = 0;
        int linker_insert = 0;
        int parsed_insert = 0;
        
        for(File inputFile: fileList)
        {
            Scanner fin = new Scanner(inputFile);
            fin.useDelimiter("\r");
            
            String id = "";                     //id of oligo
            String orientation = "";            //orientation of oligo (SENSE/ANTISENSE)
            String exp = "12th";                //experiment number
            String platform = "454";            //sequencing platform
            Oligo oligo = new Oligo();          //current read sequence
            Oligo parsedLeft = new Oligo();     //left parsed oligo
            Oligo parsedRight = new Oligo();    //right parsed oligo
            Oligo parsedBoth = new Oligo();     //left and right parsed oligo
            Oligo insert = new Oligo();         //both_parsed/left_parsed/right_parsed/unparsed oligo
            Oligo leftLinker = new Oligo();
            Oligo rightLinker = new Oligo();
            Oligo leftQuery = new Oligo();
            Oligo rightQuery = new Oligo();
            Oligo MmeISite = new Oligo("TCCGAC");
            Oligo _5primeLinker454A1_1 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATTCTCTCCGAC");     //454 5' SENSE linker
            Oligo _3primeLinker454B1_3P = new Oligo("ATCACCAGCTGCTCATAGACTGTACCAGTCTGAGCGGGCTGGCAAGGC");    //454 3' SENSE linker
            Oligo _5primeQuery454A1_1 = new  Oligo("TCAAGTCGTTCCTCATTCTCTCCGAC");   //454 5' SENSE query
            Oligo _3primeQuery454B1_3P = new Oligo("ATCACCAGCTGCTCATAGACTGTACC");   //454 3' SENSE query
            Oligo bad = new Oligo("GGGTTCCTCATTCTCTCCGAC");      //bad insert sequence to be excluded
            int min5primeKeyLength = leftQuery.length() - 2;
            int min3primeKeyLength = rightQuery.length() - 2;
            int allowedMismatches = 2;
            int ins = 1;
            int del = 1;
            int line = 3994;    //worm line used
            
            System.out.println("Parsing " + inputFile.getName() + " ...");
            System.out.println();
            
            while(fin.hasNext())
            {
                try
                {
                    id = fin.next().substring(0, 17);  //read id
                    oligo = new Oligo(fin.next()); //read corresponding oligo
                    
                    if(isPlatform(oligo, _5primeLinker454A1_1, _3primeLinker454B1_3P) && oligo.toString().contains("TCCGAC"))
                    {
                        orientation = getSense(oligo, _5primeLinker454A1_1, _3primeLinker454B1_3P);
                        
                        if(!orientation.contains("ANTI") && !orientation.contains("UNKNOWN"))
                        {
                            leftLinker = _5primeLinker454A1_1;
                            rightLinker = _3primeLinker454B1_3P;
                            leftQuery = _5primeQuery454A1_1;
                            rightQuery = _3primeQuery454B1_3P;
                            
                            if(orientation == "SENSE")
                            {
                                parsedLeft = parseLeftLinker(oligo, leftLinker, leftQuery, allowedMismatches, ins, del, min5primeKeyLength);
                                parsedBoth = parseRightLinker(parsedLeft, rightLinker, rightQuery, allowedMismatches, ins, del, min3primeKeyLength);
                                insert = parsedBoth;
                                //System.out.println("SENSE both: " + insert);
                            }
                            else if(orientation == "SENSE_5prime")
                            {
                                parsedLeft = parseLeftLinker(oligo, leftLinker, leftQuery, allowedMismatches, ins, del, min5primeKeyLength);
                                insert = parsedLeft;
                                //System.out.println("SENSE 5': " + insert);
                            }
                            else if(orientation == "SENSE_3prime")
                            {
                                parsedRight = parseRightLinker(oligo, rightLinker, rightQuery, allowedMismatches, ins, del, min3primeKeyLength);
                                insert = parsedRight;
                                //System.out.println("SENSE 3': " + insert);
                            }
                            
                            if(insert.length() > 0 && !insert.isFuzzySmithWatermanMatch(bad, 2, 1, 1) && !insert.isFuzzySmithWatermanMatch(bad.antiparallel(), 2, 1, 1))
                            {
                                parsed_insert++;
//                                parsed.println(id + '\t' + insert + '\t' + orientation + '\t' + line + '\t' + platform + '\t' + exp);
                            }
                            else if(insert.length() > 0 && (insert.isFuzzySmithWatermanMatch(bad, 2, 1, 1) || insert.isFuzzySmithWatermanMatch(bad.antiparallel(), 2, 1, 1)))
                                linker_insert++;
                        }
                        else if(orientation.contains("ANTI"))
                        {
                            oligo = oligo.reverse();
                            MmeISite = MmeISite.complement();
                            leftLinker = _5primeLinker454A1_1.complement();
                            rightLinker = _3primeLinker454B1_3P.complement();
                            leftQuery = _5primeQuery454A1_1.complement();
                            rightQuery = _3primeQuery454B1_3P.complement();
                            
                            if(orientation == "ANTISENSE")
                            {
                                parsedLeft = parseLeftLinker(oligo, leftLinker, leftQuery, allowedMismatches, ins, del, min5primeKeyLength);
                                parsedBoth = parseRightLinker(parsedLeft, rightLinker, rightQuery, allowedMismatches, ins, del, min3primeKeyLength);
                                insert = parsedBoth.reverse();
                                //System.out.println("antiSENSE both: " + insert);
                            }
                            else if(orientation == "ANTISENSE_5prime")
                            {
                                parsedLeft = parseLeftLinker(oligo, leftLinker, leftQuery, allowedMismatches, ins, del, min5primeKeyLength);
                                insert = parsedLeft.reverse();
                                //System.out.println("antiSENSE 5': " + insert);
                            }
                            else if(orientation == "ANTISENSE_3prime")
                            {
                                parsedRight = parseRightLinker(oligo, rightLinker, rightQuery, allowedMismatches, ins, del, min3primeKeyLength);
                                insert = parsedRight.reverse();
                                //System.out.println("antiSENSE 3': " + insert);
                            }
                            
                            if(insert.length() > 0 && !insert.isFuzzySmithWatermanMatch(bad, 2, 1, 1) && !insert.isFuzzySmithWatermanMatch(bad.antiparallel(), 2, 1, 1))
                            {
                                parsed_insert++;
//                                parsed.println(id + '\t' + insert + '\t' + orientation + '\t' + line + '\t' + platform + '\t' + exp);
                            }
                            else if(insert.length() > 0 && (insert.isFuzzySmithWatermanMatch(bad, 2, 1, 1) || insert.isFuzzySmithWatermanMatch(bad.antiparallel(), 2, 1, 1)))
                                linker_insert++;
                        }
                        else if(orientation == "UNKNOWN")
                        {
                            //System.out.println("UNKNOWN");
                            unknown_seq++;
                            insert = oligo;
                            orientation = "UNKNOWN";
//                            unknown.println(id + '\t' + insert + '\t' + orientation + '\t' + platform + '\t' + exp);
                        }
                    } // end if(isPlatform()...)
                    else
                    {
                        //System.out.println("not mine");
                        notMine++;
//                        not_mine.print(id + "\r");
//                        not_mine.print(insert + "\r");
                    }
                } //end try{}
                
                catch (OligoException error)
                {
                    //System.out.println("catch (OligoException error)");
                    unknown_seq++;
                    insert = oligo;
                    orientation = "N/A";
//                    unknown.println(id + '\t' + insert + '\t' + orientation + '\t' + platform + '\t' + exp + '\t' + error);
                }
            } //end while() loop
            fin.close();
        } //end for loop
//        parsed.close();
//        unknown.close();
//        not_mine.close();
        
        System.out.println("parsed inserts: " + parsed_insert);
        System.out.println("linker_inserts: " + linker_insert);
        System.out.println("unknown sequences: " + unknown_seq);
        System.out.println("not my sequences: " + notMine);
        System.out.println("Finished parsing all 454 reads.");
        
//        Oligo test = new Oligo("AGCTTAGGCTTTGATGGTGCCTACAGTTTAGTGAGGAATTCCGTTGGCTGTTTAGTGAGGAATTCCGTTGGCTGACGTTACCGTCTGATGGCGCGAGGGAGGC");
//                                 //GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATTCTCTCCGAC
//        try
//        {
//            System.out.println("foundLeftLinker: " + foundLeftLinker(test.reverse(), _5primeLinker454A1_1.complement()));
//            System.out.println("foundRightLinker: " + foundRightLinker(test, _3primeLinker454B1_3P.complement()));
//            System.out.println("isPlatForm: " + isPlatform(test, _5primeLinker454A1_1, _3primeLinker454B1_3P));
//            System.out.println("getSense: " + getSense(test, _5primeLinker454A1_1, _3primeLinker454B1_3P));
//            //System.out.println("sw: " + _5primeQuery454A1_1.smithWaterman(test, 1, 1));
//            //System.out.println("parseLeftLinker: " + parseLeftLinker(test, _5primeLinker454A1_1, _5primeQuery454A1_1, MmeISite, 3, 1, 1, 20));
//            System.out.println("parseRightLinker: " + parseRightLinker(test, _3primeLinker454B1_3P, _3primeQuery454B1_3P, 3, 1, 1, 20));
//        }
//        catch (OligoException error)
//        {
//            System.out.println(error);
//        }
    } //end main()
    
    
    
    //============================= parse454 accessory methods ==============================//
    /**
     *Assumes insert is in 5'-->3' orientation
     *Parses out the left linker if anchor insert (MmeIsite) is found. Will parse out entire left linker including the MmeI site
     */
    private static Oligo parseLeftLinker(Oligo oligo, Oligo leftLinker, Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        Oligo foundKey = new Oligo();
        
        try
        {
            foundKey = oligo.getFuzzySWkey(query, mismatches, ins, del, minKeyLength); //swKey doesn't include anchor seq'
            if(oligo.isFuzzyMatch(foundKey, 0))
                return oligo.exciseLeftFrom(foundKey, 0);
            else throw new OligoException("MatchNotFound", "parseLeftLinker()");
        }
        catch (OligoException error)
        {
            throw error;
        }
    } //end parseLeftLinker() method
    
    
    
    private static Oligo parseRightLinker(Oligo oligo, Oligo rightLinker, Oligo query, int mismatches, int ins, int del, int minKeyLength) throws OligoException
    {
        Oligo foundKey = new Oligo();
        
        try
        {
            foundKey = oligo.getFuzzySWkey(query, mismatches, ins, del, minKeyLength);
            if(oligo.isFuzzyMatch(foundKey, 0))
                return oligo.exciseRightFrom(foundKey, 0);
            else throw new OligoException("MatchNotFound", "parseRightLinker()");
        }
        catch (OligoException error)
        {throw error;}
    } //end parseRightLinker()
    
    
    
    //---------------------------------- getSense() -----------------------------------//
    private static String getSense(Oligo oligo, Oligo _5primeLinker, Oligo _3primeLinker) throws OligoException
    {
        String sense = new String();
        
        try
        {
            boolean SENSE_found5PrimeLinker = foundLeftLinker(oligo, _5primeLinker);
            boolean SENSE_found3PrimeLinker = foundRightLinker(oligo, _3primeLinker);
            boolean ANTISENSE_found5PrimeLinker = foundLeftLinker(oligo, _3primeLinker.antiparallel());
            boolean ANTISENSE_found3PrimeLinker = foundRightLinker(oligo, _5primeLinker.antiparallel());
            
            if(SENSE_found5PrimeLinker && SENSE_found3PrimeLinker)
                sense = "SENSE";
            else if (SENSE_found5PrimeLinker ^ SENSE_found3PrimeLinker)
                sense = SENSE_found5PrimeLinker ? "SENSE_5prime" : "SENSE_3prime";
            else if(ANTISENSE_found5PrimeLinker && ANTISENSE_found3PrimeLinker)
                sense = "ANTISENSE";
            else if (ANTISENSE_found5PrimeLinker ^ ANTISENSE_found3PrimeLinker)
                sense = ANTISENSE_found5PrimeLinker ? "ANTISENSE_5prime" : "ANTISENSE_3prime";
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
            boolean SENSE_found5PrimeLinker = foundLeftLinker(oligo, _5primeLinker);
            boolean SENSE_found3PrimeLinker = foundRightLinker(oligo, _3primeLinker);
            boolean ANTISENSE_found5PrimeLinker = foundLeftLinker(oligo, _3primeLinker.antiparallel());
            boolean ANTISENSE_found3PrimeLinker = foundRightLinker(oligo, _5primeLinker.antiparallel());
            
            if((SENSE_found5PrimeLinker || SENSE_found3PrimeLinker) || (ANTISENSE_found5PrimeLinker || ANTISENSE_found3PrimeLinker))
                platform = true;
            else throw new OligoException("UNKNOWN PLATFORM");
        }
        catch (OligoException error)
        {throw error;}
        
        return platform;
    } //end getSense() method
    
    
    
    //------------------------ foundLeftLinker() ----------------------------//
    private static boolean foundLeftLinker(Oligo oligo, Oligo leftLinker) throws OligoException
    {
        try
        {
            Oligo slidingQuery = new Oligo();
            final int WIDTH = 11;   //width of sliding query
            int found = 0;
            
            for(int i = 0; (i + WIDTH < leftLinker.length() - 1)  &&  found <= 3; i++)
            {
                slidingQuery = leftLinker.extractSequence(i, i + WIDTH); //scan along oligo to get key of 11 nucleotides
                if(oligo.isFuzzyMatch(slidingQuery, 1))
                    ++found;
            }
            return (found >= 3) ? true : false;
        }
        catch (OligoException error)
        {throw error;}
    } //end foundLeftLinker() method
    
    
    
    //------------------------ foundRightLinker() ----------------------------//
    private static boolean foundRightLinker(Oligo oligo, Oligo rightLinker) throws OligoException
    {
        try
        {
            Oligo slidingQuery = new Oligo();
            final int WIDTH = 11;
            int found = 0;
            
            for(int i = rightLinker.length() - 1; i > WIDTH && found <= 3; i--)
            {
                slidingQuery = rightLinker.extractSequence(i - WIDTH, i); //get key of 11 nucleotides and scan along oligo
                if(oligo.isFuzzyMatch(slidingQuery, 1))
                    ++found;
            }
            return (found >= 3) ? true : false;
        }
        catch (OligoException error)
        {throw error;}
    } //end foundLeftLinker() method
}//end El Parsero
