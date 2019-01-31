
/**
 * @author Ky Sha
 * May 14, 2010
 * Class for manipulation of half site Dam tags
 */
public class DamHS
{
	private String chr;
	private String tag;
	private String strand;
	private int gatcStart;
	private int tagStart;
	private int tagStop;


	//============================================================ CONSTRUCTORS ==========================================================//
	public DamHS()
	{
	}


	public DamHS(String hs) // minimal DamHS object: 0=tag, 1=chromosome, 2=strand, 3=gatcStart, 4=tagStart, 5=tagStop
	{
		if(!hs.contains("#"))
		{
			String[] s = hs.split("\t");
			this.tag = s[0];
			this.chr = s[1];
			this.strand = s[2];
			this.gatcStart = Integer.parseInt(s[3]);
			this.tagStart = Integer.parseInt(s[4]);
			this.tagStop = Integer.parseInt(s[5]);
		}
		else if(hs.contains("#"))
		{
			System.out.println("ERROR: [" + Thread.currentThread().getStackTrace()[1].getClassName() + "." + Thread.currentThread().getStackTrace()[1].getMethodName() + "]: Input text is a comment line. DamHS object not created. Thread aborted.");
			System.exit(0);
		}
	} //end constructor Dam()


	//================================================= METHODS: object manipulation =====================================================//
	public String getChr()
	{
		return this.chr;
	}


	public int getGatcStart()
	{
		return this.gatcStart;
	}


	public Coordinate getGatcCoordinate()
	{
		return new Coordinate(this.chr, this.gatcStart, this.gatcStart);
	}


	public String getStrand()
	{
		return this.strand;
	}


	public String getTag()
	{
		return this.tag;
	}


	public int getTagStart()
	{
		return this.tagStart;
	}


	public int getTagStop()
	{
		return this.tagStop;
	}


	public Coordinate getTagCoordinate()
	{
		return new Coordinate(this.chr, this.tagStart, this.tagStop);
	}


//	public String toBED9(String itemRgb)
//	{
//		String tagName = "(" + this.strand + ")" + this.gatcStart;
//		String useScore = "0";
//		return (this.chr + "\t" + this.tagStart + "\t" +  this.tagStop + "\t" + tagName + "\t" + useScore + "\t" + this.strand + "\t" + this.tagStart + "\t" +  this.tagStop + "\t" + itemRgb);
//	}
	public String toBED9(String tagName, String itemRgb)
	{
		return (this.chr + "\t" + this.tagStart + "\t" + this.tagStop + "\t" + tagName + "\t" + "0" + "\t" + this.strand + "\t" + this.tagStart + "\t" + this.tagStop + "\t" + itemRgb);
	}


	/**
	 * Returns the plus strand sequence
	 * @return
	 */
	@Override
	public String toString()
	{
		return this.tag + "\t" + this.chr + "\t" + this.strand + "\t" + this.gatcStart + "\t" + this.tagStart + "\t" + this.tagStop;
	}


	//======================================================= METHODS: boolean ===========================================================//
	public boolean isEmpty()
	{
		return (this.tag == null || this.chr == null || this.strand == null) ? true : false;
	}


	/**
	 * Determines if two half sites belong to the isIdenticalTo GATC site.
	 * @param hs1 half site 1
	 * @param hs2 half site 2
	 * @return boolean
	 */
	public boolean sameGatcSite(DamHS hs2)
	{
		String chr2 = hs2.getChr();
		String strand2 = hs2.getStrand();
		int gatcStart2 = hs2.getGatcStart();

		boolean sameSite = false;
		boolean sameChr = this.chr.equalsIgnoreCase(chr2);
		boolean oppositeSigns = (this.strand.equals("+") && strand2.equals("-")) || (this.strand.equals("-") && strand2.equals("+"));

		if(sameChr && oppositeSigns)
		{
			boolean cond1 = this.strand.equals("+") && (this.gatcStart - gatcStart2 == -3);
			boolean cond2 = this.strand.equals("-") && (this.gatcStart - gatcStart2 == 3);

			sameSite = (cond1 || cond2) ? true : false;
		}
		return sameSite;
	}
	//((((((((((((((((((((((((((((((((((((((((((((((((((((((( Comparators )))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))//
	/**
	 * Comparator for Dam half site, defined only by gatcStart
	 */
//	public static Comparator<DamHS> damHSComparator = new Comparator<DamHS>()
//	{
//		public int compare(DamHS hs1, DamHS hs2)
//		{
//			if(!hs1.chr().equalsIgnoreCase(hs2.chr()))
//			{
//				System.out.println("ERROR: [" + Thread.currentThread().getStackTrace()[1].getClassName() + "." + Thread.currentThread().getStackTrace()[1].getMethodName() + "]: Dam half sites are on different chromosomes. Thread aborted.");
//				System.exit(0);
//			}
//
//			Coordinate coord1 = hs1.getGatcCoordinate();
//			Coordinate coord2 = hs2.getGatcCoordinate();
//
//			int startComparator = coord1.start().compareTo(coord2.start());
//			int stopComparator = coord1.stop().compareTo(coord2.stop());
//
//			return (startComparator != 0) ? startComparator : stopComparator;
//		}
//	};
} //public class DamHS

