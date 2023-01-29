// Java program to parse the molecule
// and get the atoms count

import java.io.*;
import java.util.*;

class GFG {

	static void getCount(String str)
	{

		// Use LinkedHashmap to store
		// elements in insertion order
		Map<String, Integer> mp
			= new LinkedHashMap<String, Integer>();

		for (int i = 0; i < str.length(); i++) {
			int count = 0;

			// Convert the string element into character
			char c = str.charAt(i);
			String a = " ";

			// Convert the character element
			// into string element
			a = String.valueOf(c);

			// Check string contains the Capital A - Z value.
			if (a.matches("[A-Z]")) {

				for (int j = i + 1; j < str.length(); j++) {

					char d = str.charAt(j);
					String b = String.valueOf(d);

					// Check string contains the small a-z value.
					if (b.matches("[a-z]")) {

						a += b;

						if (mp.get(a) == null)
							mp.put(a, 1);
						else
							mp.put(a, mp.get(a) + 1);
						count = 1;
					}

					// Check string contains the number value.
					else if (b.matches("[\\d]")) {

						int k = Integer.parseInt(b);
						mp.put(a, k);
						count = 1;
					}

					else {
						i = j - 1;
						break;
					}
				}

				if (count == 0) {
					if (mp.get(a) == null)
						mp.put(a, 1);
					else
						mp.put(a, mp.get(a) + 1);
				}
			}
		}

		System.out.println("\nAtom count:");
		for (Map.Entry<String, Integer>
				entry : mp.entrySet())
			System.out.println(entry.getKey()
							+ " "
							+ entry.getValue());

		mp.clear();
	}

	// Driver code
	public static void main(String[] args)
	{
		String str = "Fe2H3OH";
		System.out.println("Given molecule: "
						+ str);
		getCount(str);
	}
}
