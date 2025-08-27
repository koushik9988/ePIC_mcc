#include "iniparser.h"


std::map<std::string, std::map<std::string, std::string>> INIParser::parse(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::map<std::string, std::map<std::string, std::string>> sections;
    std::string currentSection;

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == ';' || line[0] == '#')
        {
            continue;
        }

        if (line[0] == '[' && line[line.length() - 1] == ']')
        {
            currentSection = line.substr(1, line.length() - 2);
            continue;
        }

        if (currentSection == "species")
        {
            // Handle species section differently
            sections[currentSection][std::to_string(sections[currentSection].size())] = line;
        }
        else
        {
            std::istringstream iss(line);
            std::string key, value;
            if (std::getline(iss, key, '=') && std::getline(iss, value))
            {
                trim(key);
                trim(value);
                sections[currentSection][key] = value;
            }
        }
    }

    return sections;
}



int INIParser::getInt(const std::map<std::string, std::string> &section, const std::string &key) 
{
    return std::stoi(section.at(key));
}

double INIParser::getDouble(const std::map<std::string, std::string> &section, const std::string &key) 
{
    return std::stod(section.at(key));
}

std::string INIParser::getString(const std::map<std::string, std::string> &section, const std::string &key) 
{
    return section.at(key);
}

void INIParser::trim(std::string &str)
{
    if (str.empty()) 
    {
        return; //Return early if string is empty
    }
    
    size_t first = str.find_first_not_of(" \t");
    size_t last = str.find_last_not_of(" \t");

    if (first == std::string::npos) {
        str.clear(); // Entire string is whitespace; clear the string
    } 
    else 
    {
        str = str.substr(first, (last - first + 1));
    }
}

// Function to split a string by a delimiter
std::vector<std::string> INIParser::split(const std::string &str, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) 
    {
        tokens.push_back(token);
    }
    return tokens;
}

/*
std::tuple<std::string, int, double> INIParser::loadtypeextract(const std::string &position_init)
{
    size_t open = position_init.find('(');
    size_t close = position_init.find(')', open);

    if (open != std::string::npos && close != std::string::npos && close > open)
    {
        std::string prefix = position_init.substr(0, open);    // e.g., "1.2sin"
        std::string n_str = position_init.substr(open + 1, close - open - 1);  // e.g., "3"

        std::string func_type;
        double amplitude = 1.0;

        if (prefix.size() >= 3)
        {
            std::string last3 = prefix.substr(prefix.size() - 3);
            if (last3 == "sin" || last3 == "cos")
            {
                func_type = last3;  // "sin" or "cos"
                std::string amp_str = prefix.substr(0, prefix.size() - 3);  // e.g., "1.2"

                if (!amp_str.empty())
                {
                    amplitude = std::stod(amp_str);
                }
            }
            else
            {
                func_type = prefix;
            }
        }
        else
        {
            func_type = prefix;
        }

        int n = std::stoi(n_str);

        return {func_type, n, amplitude};
    }

    return {position_init, 0, 1.0};
}

*/



///////////////
std::tuple<std::string, int, double> INIParser::loadtypeextract(const std::string &position_init)
{
    size_t open = position_init.find('(');
    size_t close = position_init.find(')', open);

    if (open != std::string::npos && close != std::string::npos && close > open)
    {
        std::string prefix = position_init.substr(0, open);
        std::string arg_str = position_init.substr(open + 1, close - open - 1);

        if (prefix == "extend")
        {
            // Format: extend(start,end) --> we'll store start in `n`, end in `amplitude`
            size_t comma = arg_str.find(',');
            if (comma != std::string::npos)
            {
                int start = std::stoi(arg_str.substr(0, comma));
                double end = std::stod(arg_str.substr(comma + 1));
                return {"extend", start, end};  // using int for start, double for end
            }
        }

        // Normal sin, cos, etc.
        std::string func_type;
        double amplitude = 1.0;

        if (prefix.size() >= 3)
        {
            std::string last3 = prefix.substr(prefix.size() - 3);
            if (last3 == "sin" || last3 == "cos")
            {
                func_type = last3;
                std::string amp_str = prefix.substr(0, prefix.size() - 3);

                if (!amp_str.empty())
                    amplitude = std::stod(amp_str);
            }
            else
            {
                func_type = prefix;
            }
        }
        else
        {
            func_type = prefix;
        }

        int n = std::stoi(arg_str);
        return {func_type, n, amplitude};
    }

    return {position_init, 0, 1.0};
}

////////////////

std::vector<std::pair<int, int>> INIParser::parseCollGroup(const std::string &collGroupStr)
{
    std::vector<std::pair<int, int>> pairs;
    std::vector<std::string> tokens = split(collGroupStr, ',');

    /*
    ASCII values
    0 - 48
    1 - 49
    2 - 50
    3 - 51
    4 - 52
    so on
    */

    for (const std::string& token : tokens)
    {
        if (token.size() == 2)
        {  
            int first = token[0] - '0';  // Convert char to int
            int second = token[1] - '0';
            pairs.emplace_back(first, second);
        }
        else
        {
            std::cerr << "Invalid collision pair token: " << token << std::endl;
        }
    }

    return pairs;
}

