/*
 * 
 */
import java.util.*;
import java.io.*;
public class FileIntersect {
	static boolean local = false;
public static void main(String[] args) throws IOException
{
	if(!local && args.length != 2)
	{
		System.out.println("Usage: java FileIntersect <input1> <input2>");
		return;
	}
	String fn1 = "/home/mkirsche/ecoli/oxford_against_ref.paf.filtered";
	String fn2 = "/home/mkirsche/ecoli/oxford_filtered_miniasmclone.fasta";
	if(args.length > 0)
	{
		fn1 = args[0];
		fn2 = args[1];
	}
	Scanner input1 = new Scanner(new FileInputStream(new File(fn1)));
	Scanner input2 = new Scanner(new FileInputStream(new File(fn2)));
	
	HashSet<String> set1 = new HashSet<String>();
	int n1 = 0, n2 = 0;
	while(input1.hasNext())
	{
		n1++;
		set1.add(input1.nextLine().split(" \t")[0]);
	}
	int count = 0;
	while(input2.hasNext())
	{
		n2++;
		String s = input2.nextLine().split(" \t")[0];
		if(set1.contains(s)) count++;
	}
	System.out.println("File 1 has: " + n1 + " reads");
	System.out.println("File 2 has: " + n2 + " reads");
	System.out.println("Intersection: " + count);
}
}
