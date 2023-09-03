package bulkatac2fragments;

import java.io.IOException;
import java.util.StringTokenizer;

public class Fragment implements Comparable<Fragment>{
	
	String chr;
	int start;
	String line;
	
	Fragment(String line) throws IOException {
		this.line=line;
		StringTokenizer stok=new StringTokenizer(line, "\t");
		chr=stok.nextToken();
		start=Integer.parseInt(stok.nextToken());

	}

	@Override
	public int compareTo(Fragment o) {
		int comp=chr.compareTo(o.chr);
		if(comp==0) {
			return Integer.compare(start, o.start);
		} else
			return comp;
	}
	
}