



V from_str(std::string line, std::string sep = ", ")
{
	// Tokenize (Split)
	std::vector<std::string> tokens; 
	int start = 0;
	int end = line.find(sep);
	while (end != -1) 
	{
	    tokens.emplace_back(s.substr(start, end - start));
	    start = end + sep.size();
	    end = line.find(sep, start);
	}
	tokens.emplace_back(s.substr(start, end - start));


	// Filter
	std::remove_if(tokens.begin(), tokens.end(), !VUtils::is_num);
	//std::remove_if(tokens.begin(), tokens.end(), 
	//			   [](std::string l) -> bool {return !VUtils::is_num(l)});

	// Parse the first three strings to double, ignores the rest.
	V converted(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
	return converted;
}

std::vector<V> from_str_multi(std::string lines, std::string sep = ", ")
{

	//Split string by endlines
	// Tokenize (Split)
	std::vector<std::string> tokens; 
	int start = 0;
	int end = line.find(sep);
	while(end != -1) 
	{
	    tokens.emplace_back(s.substr(start, end - start));
	    start = end + sep.size();
	    end = line.find(sep, start);
	}
	tokens.emplace_back(s.substr(start, end - start));


	// Filter
	std::remove_if(tokens.begin(), tokens.end(), !VUtils::is_num);
	//std::remove_if(tokens.begin(), tokens.end(), 
	//			   [](std::string l) -> bool {return !VUtils::is_num(l)});

	
	//Calculate the total number of parsable Vs.
	unsigned int no_vecs = tokens.size() / V::cardinality;
	std::vector<V> converted;
	converted.reserve(no_vecs);

	// Keep parsing until the converted list is exhausted.
	for(unsigned int i = 0; i < no_vecs; i+=3)
	{
		converted.emplace_back(V(std::stod(tokens[i]), 
								 std::stod(tokens[i+1]), 
								 std::stod(tokens[i+2])));
	}

	return converted;


}