public class HW5 {
	public static void main(String args[])
	{
		// my tests
		int i = LongestAscendingSubsequence("182007");
		//System.out.println(i);
		
		
		String s = findpanlindrom("nixonnoxin");
		String s2 = findpanlindrom("asdfnixonnoxinfdsa");
		// code does not work for this case and not sure why
		String s3 = findpanlindrom("asdfnixonnoxinasdf");
		String s4 = findpanlindrom("nixonnoxinasdf");
		//System.out.println(s);
		//System.out.println(s2);
		//System.out.println(s3);
		//System.out.println(s4);
	}
	
	public static int LongestAscendingSubsequence(String unorderedInput)
	{
		//implement mergeSort
		String orderedInput = MergeSort(unorderedInput);
		//System.out.println(orderedInput);
		return LongestCommonSubsequence(unorderedInput, orderedInput);	
	}
	public static String MergeSort(String unorderedInput)
	{
		//add your code here
		// found some mergeSort algorithm and code and analyzed it. Then converted it to something
		// I could use for this assignment.
		// http://javahungry.blogspot.com/2013/06/java-sorting-program-code-merge-sort.html
		String merged_string = "";
		String merged_string1 = "";
		String merged_string2 = "";
		
		// converts string into array of ints
		int[] unorderedInts = new int[unorderedInput.length()/2];
		//System.out.println(unorderedInts.length);
		for (int i = 0; i < unorderedInts.length; i++)
		{
			//System.out.println(i);
			unorderedInts[i] = (unorderedInput.charAt(i * 2) - 48) * 10;
			unorderedInts[i] = unorderedInts[i] + (unorderedInput.charAt((i * 2) + 1) - 48);
			//System.out.println(unorderedInts[i]);
		}
		//System.out.println(unorderedInts.toString());
		// start of mergeSort algorithm
		if (unorderedInts.length <= 1)
		{
			merged_string = "" + unorderedInts[0];
			if (merged_string.length() == 1)
			{
				merged_string = 0 + merged_string;
			}
			return merged_string;
		}
		
		// Split the array in half
		int[] first = new int[unorderedInts.length/2];
		int[] second = new int[unorderedInts.length - first.length];
		System.arraycopy(unorderedInts, 0, first, 0, first.length);
		System.arraycopy(unorderedInts, first.length, second, 0, second.length);
//		for(int i = 0; i < second.length; i++)
//		{
//			System.out.println(second[i]);
//		}
		
		// sort each half
		String first_string = "";
		String second_string = "";
		// convert each array into strings
		for(int i = 0; i < first.length; i++)
		{
			String copy = "" + first[i];
			if (copy.length() == 1)
			{
				copy = 0 + copy;
			}
			first_string += copy;
		}
		//System.out.println(first_string);
		for(int j = 0; j < second.length; j++)
		{
			String copy = "" + second[j];
			if (copy.length() == 1)
			{
				copy = 0 + copy;
			}
			second_string += copy;
		}
		//System.out.println(second_string);
		// send to MergeSort method to merge them and then convert to array
		merged_string1 = MergeSort(first_string);
		for (int i = 0; i < first.length; i++)
		{
			//System.out.println(i);
			first[i] = (merged_string1.charAt(i * 2) - 48) * 10;
			first[i] = first[i] + (merged_string1.charAt((i * 2) + 1) - 48);
			//System.out.println(unorderedInts[i]);
		}
		merged_string2 = MergeSort(second_string);
		for (int i = 0; i < second.length; i++)
		{
			//System.out.println(i);
			second[i] = (merged_string2.charAt(i * 2) - 48) * 10;
			second[i] = second[i] + (merged_string2.charAt((i * 2) + 1) - 48);
			//System.out.println(unorderedInts[i]);
		} 
		
		// Merge the halves together, overwriting the original array
		int[] mergedArray = merge(first, second);
		for(int i = 0; i < mergedArray.length; i++)
		{
			String copy = "" + mergedArray[i];
			if (copy.length() == 1)
			{
				copy = 0 + copy;
			}
			merged_string += copy;
		}
		return merged_string;
	}
	public static int LongestCommonSubsequence(String S1, String S2)
	{
		//TrackList is used to track the subsequence. 
		//TrackList[i][j]=1 represents F[i][j] comes from F[i-1][j-1].
		//TrackList[i][j]=2 represents F[i][j] comes from F[i-1][j].
		//TrackList[i][j]=3 represents F[i][j] comes from F[i][j-1].
		int[][] TrackList = new int[S1.length() / 2 + 1][S2.length() / 2 + 1];
		for(int i = 0; i <= S1.length() / 2; i++)
			for(int j = 0; j <= S2.length() / 2; j++)
				TrackList[i][j] = 0;
		
		//F is used to store the intermediate results.
		// I used to store values to use for backtracking and printing sequence
		int[][] F = new int[S1.length() / 2 + 1][S2.length() / 2 + 1];
		for(int i = 0; i <= S1.length() / 2; i++)
			for(int j = 0; j <= S2.length() / 2; j++)
				F[i][j] = 0;
		
		// Bridge contents : These are the values that are equal and make a bridge in the table
		int[][] Bridge_content = new int[S1.length() / 2 + 1][S2.length() / 2 + 1];
		for(int i = 0; i <= S1.length() / 2; i++)
			for(int j = 0; j <= S2.length() / 2; j++)
				Bridge_content[i][j] = 0;
	
		//add your code here.
		// I found an algorithm for LCS online off a youtube video and implemented that
		// https://www.youtube.com/watch?v=P-mMvhfJhu8
		//System.out.println(S1);
		//System.out.println(S2);
		// builds table for longest subsequence. number for longest subsequence is located
		// in last array at bottom right corner
		for(int i = 1; i <= S1.length() / 2; i++){
			for(int j = 1; j <= S2.length() / 2; j++) {
				if (S1.charAt(2 * j - 2) == S2.charAt(2 * i - 2) && S1.charAt(2 * j - 1) == S2.charAt(2 * i - 1)){
					TrackList[i][j] = TrackList[i - 1][j - 1] + 1;
					F[i][j] = F[i - 1][j - 1] + 1;
					//System.out.print(i);
					//System.out.println(S1.charAt(i));
					//System.out.println(j);
					//System.out.println(S2.charAt(j));
				}
				else {
					TrackList[i][j] = Math.max(TrackList[i][j - 1], TrackList[i - 1][j]);
					F[i][j] = Math.max(F[i][j - 1], F[i - 1][j]);
				}
			}
		}
		
		// updates backtracking array with bridges at specific locations
		for(int i = 1; i <= S1.length() / 2; i++){
			for(int j = 1; j <= S2.length() / 2; j++) {
				if (S1.charAt(2 * j - 2) == S2.charAt(2 * i - 2) && S1.charAt(2 * j - 1) == S2.charAt(2 * i - 1)){
					F[i][j] = -1;
					Bridge_content[i][j] = (S1.charAt(2 * j - 2) - 48) * 10;
					Bridge_content[i][j] = Bridge_content[i][j] + (S1.charAt(2 * j - 1) - 48);
				}
			}
		}
		
//		for(int i = 0; i <= S1.length() / 2; i++){
//			for(int j = 0; j <= S2.length() / 2; j++) {
//				//System.out.print(i);
//				//System.out.print(j);
//				System.out.print(" " + Bridge_content[i][j] + " ");
//			}
//			System.out.println();
//		}
		//return the length of the subsequence and print out the longest common subsequence.
		int m = F[0].length - 1;
		int n = F.length - 1;
		String subsequence = "";
		while(m != 0 && n != 0)
		{
			if (F[m][n] == -1)
			{
				subsequence = Bridge_content[m][n] + subsequence;
				m = m - 1;
				n = n - 1;
			} else {
				if (TrackList[m - 1][n] == TrackList[m][n])
				{
					m = m - 1;
				}
				if (TrackList[m][n - 1] == TrackList[m][n])
				{
					n = n - 1;
				}
			}
		}
		System.out.println(subsequence);
		return TrackList[S1.length() / 2][S2.length() / 2];
	}
	
	public static String findpanlindrom(String s) {
		  // Found online an algorithm to find the longest palindromic subsequence and
		  // implemented it
		  // http://leetcode.com/2011/11/longest-palindromic-substring-part-i.html
		  int n = s.length();
		  int longestPanBegin = 0;
		  int maxPanLen = 1;
		  // used to determine the characters that are true for palindrom (if a palindrom
		  // does exist then its symmetrical on the diagonal)
		  boolean[][] table = new boolean[1000][1000];
		  
		  // initialize diagonal to all true
		  for (int i = 0; i < n; i++) {
		    table[i][i] = true;
		  }
		  
		  // find character side by side that is same.
		  for (int j = 0; j < n-1; j++) {
		    if (s.charAt(j) == s.charAt(j+1)) {
		      table[j][j+1] = true;
		      longestPanBegin = j;
		      maxPanLen = 2;
		    }
		  }
		  
		  // find the longest palindrom that comes out of that given start of palindrom
		  for (int len = 3; len <= n; len++) {
		    for (int k = 0; k < n-len+1; k++) {
		      int l = k+len-1;
		      if (s.charAt(k) == s.charAt(l) && table[k+1][l-1]) {
		        table[k][l] = true;
		        longestPanBegin = k;
		        maxPanLen = len;
		      }
		    }
		  }
		  return s.substring(longestPanBegin, maxPanLen);
		}
	
	private static int[] merge(int[] first, int[] second)
	{
		int next_first = 0;
		int next_second = 0;
		int j = 0;
		int[] result = new int[first.length + second.length];
		
		while(next_first < first.length && next_second < second.length)
		{
			if (first[next_first] < second[next_second])
			{
				result[j] = first[next_first];
				next_first++;
			}
			else
			{
				result[j] = second[next_second];
				next_second++;
			}
			j++;
		}
		
		// copy what is left
		System.arraycopy(first, next_first, result, j, first.length - next_first);
		System.arraycopy(second, next_second, result, j, second.length - next_second);
		
		//System.out.println(first[0]);
		//System.out.println(second[0]);
		//System.out.println(result[0] +""+ result[1]);
		return result;
	}
}