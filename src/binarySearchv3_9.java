/******************************************************************************/
//binarySearch()
//Author: Ky Sha
//version 3.9
//comments: v3.9 works with all three insertion and both deletion oligos!! but
//fails to extract complete sequence from insertion4
//can assemble key fragments into complete key
/******************************************************************************/
import java.util.ArrayList;
import java.util.regex.*;
public class binarySearchv3_9
{
    public static void main(String[] args)
    {
//        Oligo target = new Oligo                                 ("GTTCCTCATTCTCTCCGAC");
//        Oligo source = new Oligo    ("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTaTCCTCATTCTCTgCCGACgatcgatcgatc");
//        Oligo insertion1 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCAaTTCTCTCCGACGAATACGTAGC");
//        Oligo insertion2 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTaTCCTCATTCTCTCCGACGAATACGTAGC");
//        Oligo insertion3 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTaTCCTCATTCTCgTCCGACGAATACGTAGC");
//        Oligo insertion4 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATaTgCTCTCCGACGAATACGTAGC");
//        Oligo insertion5 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGaTTCCTCATTCTCTCCGACGAATACGTAGC");
//        Oligo insertion6 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATTCTCTCCGAtCGAATACGTAGC");
//        Oligo deletion1 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCTCATTCTCTCCGACGAATACGTAGC");
//        Oligo deletion2 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATTCTCTCGACGAATACGTAGC");
//        Oligo deletion3 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCTTCCTCATTCTCTCCGACGAATACGTAGC");
//        Oligo wt = new Oligo       ("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATTCTCTCCGACGAATACGTAGC");
//        Oligo deletion4 = new Oligo("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCCTCATTCTCTCCGAGAATACGTAGC");
//        Oligo indel = new Oligo    ("GCCTCCCTCGCGCCATCAGAGTTCAAGTCGTTCaCTCATTCTCTCGACGAATACGTAGC");
//        Oligo temp = new Oligo("gcccccgttcaattctcccgaaagcgtccagcttcgttacccgactagctaaccggattt");
        Oligo keyA = new Oligo("GTTCAAGTCGTTCCTCATTCTC");
        Oligo F12 = new Oligo("ACGAATNGCCCTTGGCCTCCCTNGCGCCATCAGAATTCAAGTCGTTTCCTCATTCTCTCCGACTCATAATTTAAACTGTAGGAATCACCAGCTGCTCATAGACTGTACCAGTCTGAGCGGGCTGGCAAGGCAAGGGCGAATTCGCGGCCGCTAAATTCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGGTCGTTTTACAA");
        
        Oligo foundKey1 = F12.binarySearch(keyA);
        Oligo foundKey2 = F12.binarySearch(foundKey1);
        System.out.println("main(): source sequence: " + F12);
        
        
//        if(foundKey.length() < keyA.length()/2)
//            System.out.println("main(): SEARCH FAILED");
//        else    System.out.println("main(): SEARCH SUCCEEDED");
        
        System.out.println("assembled key1: " + foundKey1);
        System.out.println("assembled key2: " + foundKey2);
        System.out.println();
    } //end main()
} //end binarySearch class



