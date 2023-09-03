package bulkatac2fragments;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

public class PrioritizedReader implements Comparable<PrioritizedReader>{
	
	String chr;
	int start;
	String line;
	
	
	BufferedReader br;
	
	PrioritizedReader(File fin) throws IOException {
		br=new BufferedReader(new FileReader(fin));
		
	}
	
	boolean readLine() throws IOException {
		line=br.readLine();
		if(line!=null) {
			StringTokenizer stok=new StringTokenizer(line, "\t");
			chr=stok.nextToken();
			start=Integer.parseInt(stok.nextToken());
			return true;
		} else
			return false;
	}

	@Override
	public int compareTo(PrioritizedReader o) {
		int comp=chr.compareTo(o.chr);
		if(comp==0) {
			return Integer.compare(start, o.start);
		} else
			return comp;
	}
	
}