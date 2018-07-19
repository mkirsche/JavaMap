import java.util.*;
import java.io.*;
public class FileIntersect {
public static void main(String[] args) throws IOException
{
	String fn1 = "/home/mkirsche/ecoli/oxford_against_ref.paf.filtered";
	String fn2 = "/home/mkirsche/ecoli/oxford_filtered_miniasmclone.fasta";
	Scanner input1 = new Scanner(new FileInputStream(new File(fn1)));
	Scanner input2 = new Scanner(new FileInputStream(new File(fn2)));
	
	HashSet<String> set1 = new HashSet<String>();
	int n1 = 0, n2 = 0;
	while(input1.hasNext())
	{
		n1++;
		set1.add(input1.next());
		input1.nextLine();
	}
	int count = 0;
	while(input2.hasNext())
	{
		n2++;
		String s = input2.nextLine();
		//input2.nextLine();
		if(set1.contains(s)) count++;
	}
	System.out.println(n1+" "+n2);
	System.out.println(count);
}
}
