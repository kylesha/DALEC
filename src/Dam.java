/**
 * Class for manipulatio of FULL Dam sites
 * @author Ky Sha
 */
public class Dam extends Coordinate
{
	private String chr;
	private DamHS plusHalfSite;
	private DamHS minusHalfSite;

	//============================================================ CONSTRUCTORS ==========================================================//

	public Dam()
	{
	} //end constructor Dam()


	/**
	 * For now, Dam object MUST contain two non-empty half sites
	 * @param hs1
	 * @param hs2
	 */
	public Dam(DamHS hs1, DamHS hs2)
	{
		this.chr = hs1.getChr();

		if(hs1.sameGatcSite(hs2))
		{
			if(hs1.getStrand().equals("+"))
			{
				this.plusHalfSite = hs1;
				this.minusHalfSite = hs2;
			}
			else
			{
				this.plusHalfSite = hs2;
				this.minusHalfSite = hs1;
			}
		}
		else
		{
			System.out.println("ERROR: [" + Thread.currentThread().getStackTrace()[1].getClassName() + "." + Thread.currentThread().getStackTrace()[1].getMethodName() + "]: Input half sites [" + hs1 + "|" + hs2 + "]do not belong to the same GATC site. Thread aborted.");
			System.exit(0);
		}
	} //end constructor


	//================================================= METHODS: object manipulation =====================================================//
	@Override
	public String chr()
	{
		return this.chr;
	}


	/**
	 * Returns the Coordinate of the full site (whose start/stop are defined by the gatcStarts of plus/minus strands, respectively)
	 * @return
	 */
	public Coordinate getCoordinate()
	{
		return new Coordinate(this.chr, this.plusHalfSite.getGatcStart(), this.minusHalfSite.getGatcStart());
	}


	public DamHS getPlusHalfSite()
	{
		return this.plusHalfSite;
	}


	public DamHS getMinusHalfSite()
	{
		return this.minusHalfSite;
	}


	@Override
	public String toString()
	{
		return plusHalfSite + "|" + minusHalfSite;
	}


	//================================================== METHODS: object manipulation =====================================================//
	public boolean isProximal50(Dam dam2)
	{
		Dam dam1 = new Dam(this.plusHalfSite, this.minusHalfSite);
		Coordinate coord1 = dam1.getCoordinate();
		Coordinate coord2 = dam2.getCoordinate();
		return (coord1.interdistance(coord2) <= 50) ? true : false;
	}

		public boolean isProximal20(Dam dam2)
	{
		Dam dam1 = new Dam(this.plusHalfSite, this.minusHalfSite);
		Coordinate coord1 = dam1.getCoordinate();
		Coordinate coord2 = dam2.getCoordinate();
		return (coord1.interdistance(coord2) <= 20) ? true : false;
	}
	
	//((((((((((((((((((((((((((((((((((((((((((((((((((((((( Comparators )))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))//
	/**
	 * Orders Dam full site Coordinate objects by startCoordinate, then by stopCoordinate
	 */
//	public static Comparator<Dam> damComparator = new Comparator<Dam>()
//	{
//		public int compare(Dam dam1, Dam dam2)
//		{
//			if(!dam1.chr().equalsIgnoreCase(dam2.chr()))
//			{
//				System.out.println("ERROR: [" + Thread.currentThread().getStackTrace()[1].getClassName() + "." + Thread.currentThread().getStackTrace()[1].getMethodName() + "]: Dam objects are on different chromosomes. Thread aborted.");
//				System.exit(0);
//			}
//			Coordinate coord1 = dam1.getCoordinate();
//			Coordinate coord2 = dam2.getCoordinate();
//
//			int startComparator = coord1.start().compareTo(coord2.start());
//			int stopComparator = coord1.stop().compareTo(coord2.stop());
//
//			return (startComparator != 0) ? startComparator : stopComparator;
//		}
//	};
} //end class Dam

